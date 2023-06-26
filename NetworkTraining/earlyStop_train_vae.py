from argparse import ArgumentParser
from mmd_vae_biased_3HL import MESA_VAE
from pytorch_lightning import Trainer
from pytorch_lightning.callbacks import EarlyStopping

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)


def main(args):
    """ Main training routine specific for this project. """
    # ------------------------
    # 1 INIT LIGHTNING MODEL
    # ------------------------
    model = MESA_VAE(**vars(args))

    # ------------------------
    # 2 INIT TRAINER
    # ------------------------
    early_stopping = EarlyStopping('test_val_recon', patience=10)
    trainer = Trainer.from_argparse_args(args, callbacks=[early_stopping])

    # ------------------------
    # 3 START TRAINING
    # ------------------------
    trainer.fit(model)


def run_cli():

    parent_parser = ArgumentParser(add_help=False)
    parser = MESA_VAE.add_model_specific_args(parent_parser)

    parser = Trainer.add_argparse_args(parser)

    parser.set_defaults(gpus=1)
    parser.set_defaults(weights_summary="full")
    parser.set_defaults(max_epochs=int(5e3))
    parser.set_defaults(check_val_every_n_epoch=1)

    args = parser.parse_args()

    # ---------------------
    # RUN TRAINING
    # ---------------------
    main(args)


if __name__ == "__main__":
    run_cli()
