import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import optim
from torch.utils.data import DataLoader

from pytorch_lightning.core import LightningModule

from argparse import ArgumentParser

# Notes to myself on variables that I'm using in this code
# feature_dim = input data dimensions
# latent_dim = latent dimension
# learning_rate = network learning rate
# batch_size = 
# 


class MESA_VAE(LightningModule):
    """
    Sample model to show how to define a template.
    Example:
        >>> # define simple Net for MNIST dataset
        >>> params = dict(
        ...     in_features=28 * 28,
        ...     hidden_dim=1000,
        ...     out_features=10,
        ...     drop_prob=0.2,
        ...     learning_rate=0.001 * 8,
        ...     batch_size=2,
        ...     data_root='./datasets',
        ...     num_workers=4,
        ... )
        >>> model = LightningTemplateModel(**params)
    """

    def __init__(
        self,
        latent_dim: int,
        learning_rate: float,
        batch_size: int,
        feature_dim: int,
        hidden_dim1_encode: int,
        hidden_dim2_encode: int,
        hidden_dim3_encode: int,
        hidden_dim1_decode: int,
        hidden_dim2_decode: int,
        hidden_dim3_decode: int,
        train_percent: int,
        num_workers: int,
        epochs: int,
        save_samples_every_n_epochs: int,
        **kwargs,
    ):
        # init superclass
        super().__init__()
        # save all variables in __init__ signature to self.hparams
        self.save_hyperparameters()

        
        
        # Define encoder architecture and then call it later
        encoder_blocks = [
            nn.Linear(
                in_features=self.hparams.feature_dim,
                out_features=self.hparams.hidden_dim1_encode,
            ),
            nn.Tanh(),
            nn.BatchNorm1d(num_features=self.hparams.hidden_dim1_encode),
            nn.Linear(
                in_features=self.hparams.hidden_dim1_encode,
                out_features=self.hparams.hidden_dim2_encode,
            ),
            nn.Tanh(),
            nn.BatchNorm1d(num_features=self.hparams.hidden_dim2_encode),
            nn.Linear(
                in_features=self.hparams.hidden_dim2_encode,
                out_features=self.hparams.hidden_dim3_encode,
            ),
            nn.Tanh(),
            nn.BatchNorm1d(num_features=self.hparams.hidden_dim3_encode),
            nn.Linear(
                in_features=self.hparams.hidden_dim3_encode,
                out_features=self.hparams.latent_dim,
            ),
            nn.Tanh(),
            nn.BatchNorm1d(num_features=self.hparams.latent_dim),
        ]

        # Here is the sequential call of the the encoder architecture
        self.encoder = nn.Sequential(*encoder_blocks)

        decoder_blocks = [
            nn.Linear(
                in_features=self.hparams.latent_dim,
                out_features=self.hparams.hidden_dim1_decode
            ),
            nn.Tanh(),
            nn.BatchNorm1d(num_features=self.hparams.hidden_dim1_decode),
            nn.Linear(
                in_features=self.hparams.hidden_dim1_decode,
                out_features=self.hparams.hidden_dim2_decode
            ),
            nn.Tanh(),
            nn.BatchNorm1d(num_features=self.hparams.hidden_dim2_decode),
            nn.Linear(
                in_features=self.hparams.hidden_dim2_decode,
                out_features=self.hparams.hidden_dim3_decode
            ),
            nn.Tanh(),
            nn.BatchNorm1d(num_features=self.hparams.hidden_dim3_decode),
            nn.Linear(
                in_features=self.hparams.hidden_dim3_decode,
                out_features=self.hparams.feature_dim
            ),
            nn.Sigmoid(),
        ]

        self.decoder = nn.Sequential(*decoder_blocks)

    def encode(self, x):

        z = self.encoder(x.float())

        return (z)

    def decode(self, z):

        result = self.decoder(z)

        return result

    def forward(self, x):
        """
        No special modification required for Lightning, define it as you normally would
        in the `nn.Module` in vanilla PyTorch.
        """
        z = self.encode(x.float())

        return (z, self.decode(z))

    def compute_kernel(self, x, y):
        x_size = x.shape[0]
        y_size = y.shape[0]
        dim = x.shape[1]

        tiled_x = x.view(x_size,1,dim).repeat(1, y_size,1)
        tiled_y = y.view(1,y_size,dim).repeat(x_size, 1,1)

        return torch.exp(-torch.mean((tiled_x - tiled_y)**2,dim=2)/dim*1.0)

    def compute_mmd(self, x, y):
        x_kernel = self.compute_kernel(x, x)
        y_kernel = self.compute_kernel(y, y)
        xy_kernel = self.compute_kernel(x, y)
        return torch.mean(x_kernel) + torch.mean(y_kernel) - 2*torch.mean(xy_kernel)

    def weighted_mse_loss(self, x_pred, x_true):
        weights = torch.ones(self.hparams.feature_dim)
        weights[1:] = 1/(self.hparams.feature_dim - 1)
        batch_weights = weights.repeat(x_true.size(dim=0),1)
        batch_weights = batch_weights.to(device=torch.device("cuda:0"))
        wmse_loss = batch_weights*(x_pred - x_true)**2
        return wmse_loss.sum()/x_true.size(dim=0)

    def training_step(self, batch, batch_idx):
        """
        Lightning calls this inside the training loop with the data from the training dataloader
        passed in as `batch`.
        """
        # forward pass
        x = batch
        z, x_pred = self.forward(x.float())

        true_samples = torch.randn((len(z), self.hparams.latent_dim))
        true_samples = true_samples.to(device=torch.device("cuda:0"))

        recon_loss = F.mse_loss(x_pred, x.float())
        #recon_loss = self.weighted_mse_loss(x_pred, x.float())

        mmd_loss = self.compute_mmd(true_samples, z)

        loss = recon_loss + mmd_loss

        self.logger.experiment.add_scalars(
            "loss", {"train": loss}, global_step=self.global_step,
        )
        self.logger.experiment.add_scalars(
            "recon", {"train": recon_loss}, global_step=self.global_step,
        )
        self.logger.experiment.add_scalars(
            "MMD", {"train": mmd_loss}, global_step=self.global_step,
        )

        return {"loss": loss, "MMD": mmd_loss}

    def validation_step(self, batch, batch_idx):
        """
        Lightning calls this inside the validation loop with the data from the validation dataloader
        passed in as `batch`.
        """
        x = batch
        z, x_pred = self.forward(x.float())

        true_samples = torch.randn((len(z), self.hparams.latent_dim))
        true_samples = true_samples.to(device=torch.device("cuda:0"))

        recon_loss = F.mse_loss(x_pred, x.float())
        #recon_loss = self.weighted_mse_loss(x_pred, x.float())

        mmd_loss = self.compute_mmd(true_samples, z)

        loss = recon_loss + mmd_loss

        return {
            "val_loss": loss,
            "val_recon": recon_loss,
            "val_MMD": mmd_loss,
        }

    def validation_epoch_end(self, outputs):
        """
        Called at the end of validation to aggregate outputs.
        :param outputs: list of individual outputs of each validation step.
        """
        avg_loss = torch.stack([x["val_loss"] for x in outputs]).mean()
        avg_recon = torch.stack([x["val_recon"] for x in outputs]).mean()
        avg_MMD = torch.stack([x["val_MMD"] for x in outputs]).mean()
        self.logger.experiment.add_scalars(
            "loss", {"val": avg_loss}, global_step=self.global_step,
        )
        self.logger.experiment.add_scalars(
            "recon", {"val": avg_recon}, global_step=self.global_step,
        )
        self.logger.experiment.add_scalars(
            "MMD", {"val": avg_MMD}, global_step=self.global_step,
        )
        self.log("test_val_recon", avg_recon)

    # ---------------------
    # TRAINING SETUP
    # ---------------------
    def configure_optimizers(self):
        """
        Return whatever optimizers and learning rate schedulers you want here.
        At least one optimizer is required.
        """
        optimizer = optim.Adam(
            self.parameters(),
            lr=self.hparams.learning_rate,
        )
        # scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=10)
        return [optimizer]  # , [scheduler]

    def on_epoch_end(self):
        if self.current_epoch % self.hparams.save_samples_every_n_epochs == 0:
            print(f"Saving samples for epoch #{self.current_epoch}...")
            with torch.no_grad():
                dataiter = iter(self.val_dataloader())
                data = dataiter.next()
                batch = self.transfer_batch_to_device(data, self.device, dataiter)

                x = batch
                z, x_pred = self.forward(x.float())
        if (self.current_epoch + 1) % 1 == 0:
            print(f"Mapping latent space...")
            with torch.no_grad():
                path =''
                full_dataset = torch.from_numpy(np.load(path))
                full_dataset = full_dataset.to(device=torch.device("cuda:0"))
                z_full, x_pred_full = self.forward(full_dataset.float())
                z_full = z_full.cpu()
                np.save('LatentSpace.npy', z_full)
                x_pred_full = x_pred_full.cpu()
                np.save('ReconstructedPIVs.npy', x_pred_full)

    def prepare_data(self):
        pass

    def load_dynamics_data(self, train_percent=None):

        if train_percent is None:
            train_percent = self.hparams.train_percent

        print("Loading all data...")
        path =''
        full_dataset = torch.from_numpy(np.load(path))

        train_size = int(train_percent * len(full_dataset))
        val_size = len(full_dataset) - train_size
        print(f"Training on {train_size} examples, validating on {val_size} examples")
        ds_train, ds_val = torch.utils.data.random_split(
            full_dataset, [train_size, val_size]
        )
        return ds_train, ds_val

    def setup(self, stage):
        self.ds_train, self.ds_val = self.load_dynamics_data()

    def train_dataloader(self):

        return DataLoader(
            self.ds_train,
            batch_size=self.hparams.batch_size,
            shuffle=True,
            pin_memory=True,
            num_workers=self.hparams.num_workers,
        )

    def val_dataloader(self):
        return DataLoader(
            self.ds_val,
            batch_size=self.hparams.batch_size,
            num_workers=self.hparams.num_workers,
            pin_memory=True,
        )

    @staticmethod
    def add_model_specific_args(parent_parser):  # pragma: no-cover
        """
        Define parameters that only apply to this model
        """

        parser = ArgumentParser(parents=[parent_parser])

        # param overwrites
        # parser.set_defaults(gradient_clip_val=5.0)

        # network params
        parser.add_argument("--learning_rate", default=1e-4, type=float)
        parser.add_argument("--batch_size", default=100, type=int)
        parser.add_argument("--feature_dim", default=, type=int)
        parser.add_argument("--hidden_dim1_encode", default=64, type=int)
        parser.add_argument("--hidden_dim2_encode", default=32, type=int)
        parser.add_argument("--hidden_dim3_encode", default=16, type=int)
        parser.add_argument("--hidden_dim1_decode", default=16, type=int)
        parser.add_argument("--hidden_dim2_decode", default=32, type=int)
        parser.add_argument("--hidden_dim3_decode", default=64, type=int)
        parser.add_argument("--latent_dim", default=3, type=int)
        parser.add_argument("--train_percent", default=0.80, type=float)
        parser.add_argument("--create_dataset", default=False, type=bool)
        parser.add_argument("--num_workers", default=2, type=int)
        parser.add_argument("--epochs", default=5000, type=int)
        parser.add_argument("--save_samples_every_n_epochs", default=10000, type=int)
        return parser
