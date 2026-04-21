import torch
from torch.utils.data import DataLoader
import pandas as pd
import time
import os

# directly import from latticevision library in directory
from latticevision.device import set_device
from latticevision.img2img.dataset import (
	make_dataset,
	DataConfig,
)
from latticevision.img2img import TransUNet, UNet, ViT
from latticevision.img2img.base import (
	ModelConfig,
	NullPosEmbed,
	RotaryPosEmbed,
)
from latticevision.img2img.train import train_model, TrainingConfig
from latticevision.img2img.eval import eval_model
from latticevision.seed import set_all_random_seeds


# make sure folders exist
for d in ("results/model_wghts", "results/metrics"):
	os.makedirs(d, exist_ok=True)

# set the random seed for reproducibility
set_all_random_seeds(777)

# set the device you are on
device = set_device(
	machine="remote", gpu=True, gpu_id="cuda:0", verbose=True
)  # for remote use


modeltype_list = ["TransUNet"]
num_replicates_list = [1, 30]
pos_embeddings_all = [RotaryPosEmbed]


for num_replicates in num_replicates_list:
	# create dataset and dataloaders
	dataset_path = "data/I2I_data.h5"
	val_size = 0.1
	test_size = 0.2

	# dataset_path_i2i = "data/I2I_sample_data.h5"
	# val_size_i2i = 0.4
	# test_size_i2i = 0.5

	data_config = DataConfig(
		file_path=dataset_path,
		n_rows=128,
		n_cols=128,
		n_replicates=num_replicates,
		n_params=3,
		log_kappa2=True,
		shift_theta=True,
		val_size=val_size,
		test_size=test_size,
		random_state=777,
		verbose=True,
	)

	# load the data and time it
	start_time = time.time()
	data_dict = make_dataset(
		config=data_config,
	)
	print("--- %s seconds ---" % (time.time() - start_time))
	time_data = time.time() - start_time
	print("Dataset loading took: ", time_data / 60, " minutes.")

	train_df = data_dict["train_df"]
	val_df = data_dict["val_df"]
	test_df = data_dict["test_df"]

	n_batch = 64
	# n_batch = 2

	train_loader = DataLoader(train_df, batch_size=n_batch, shuffle=True)
	val_loader = DataLoader(val_df, batch_size=n_batch, shuffle=False)
	test_loader = DataLoader(test_df, batch_size=n_batch, shuffle=False)

	for fields, params in train_loader:
		print("Train Loader:")
		print(
			"Fields batch shape: ",
			fields.shape,
			"\nParams batch shape: ",
			params.shape,
		)
		break

	for fields, params in val_loader:
		print("Val Loader:")
		print(
			"Fields batch shape: ",
			fields.shape,
			"\nParams batch shape: ",
			params.shape,
		)
		break

	for fields, params in test_loader:
		print("Test Loader:")
		print(
			"Fields batch shape: ",
			fields.shape,
			"\nParams batch shape: ",
			params.shape,
		)
		break

	for modeltype in modeltype_list:
		# don't make unet iterate over all these pos embeddings
		if modeltype == "UNet":
			pos_embeddings_list = [NullPosEmbed]
		else:
			pos_embeddings_list = pos_embeddings_all

		for pos_embeddings in pos_embeddings_list:
			print("Running modeltype: ", modeltype)
			print("Running num_replicates: ", num_replicates)
			print("Running pos_embeddings: ", pos_embeddings.__name__)

			if modeltype == "UNet":
				config = ModelConfig(
					in_channels=num_replicates,
				)
				model = UNet(config)
			elif modeltype == "TransUNet":
				config = ModelConfig(
					in_channels=num_replicates,
					patch_size_h=2,
					patch_size_w=2,
					pos_embed_cls=pos_embeddings,
				)
				model = TransUNet(config)
			elif modeltype == "ViT":
				config = ModelConfig(
					in_channels=num_replicates,
					patch_size_h=16,
					patch_size_w=16,
					pos_embed_cls=pos_embeddings,
				)
				model = ViT(config)
			else:
				raise ValueError(
					"modeltype must be one of 'UNet', 'TransUNet', or 'ViT'"
				)

			# send model to device and count params
			model = model.to(device)
			total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
			print(f"Total number of trainable parameters: {total_params}")

			meta = {
				"model": modeltype,
				"reps": num_replicates,
				"pos": pos_embeddings.__name__,
			}

			basename = "_".join(f"{k}{v}" for k, v in meta.items())

			# train the model
			start_time = time.time()

			train_config = TrainingConfig(
				model=model,
				device=device,
				train_loader=train_loader,
				val_loader=val_loader,
				train_df=train_df,
				val_df=val_df,
				lr=1e-4,
				n_epochs=200,
				stop_patience=10,
				scheduler_patience=5,
				scheduler_factor=0.5,
				augmentation=True,
				save=True,
				save_directory="results/model_wghts/",
				savename=f"{basename}.pth",
				verbose=True,
				normalize=True,
				shuffle=True,
			)

			training_results = train_model(config=train_config)

			print("--- %s seconds ---" % (time.time() - start_time))
			time_train = time.time() - start_time
			print("Training took: ", time_train / 60, " minutes.")

			# save training results
			model = training_results["model"]
			train_losses = training_results["train_losses"]
			baseline_losses = training_results["baseline_losses"]
			val_losses = training_results["val_losses"]

			metrics_start_time = time.time()
			# evaluate the model performance
			metrics = eval_model(
				model=model,
				config=data_config,
				device=device,
				test_loader=test_loader,
				test_df=test_df,
				plot=False,
				augmentation=False,
				show=False,
			)
			print("--- %s seconds ---" % (time.time() - metrics_start_time))
			time_eval = time.time() - metrics_start_time
			print("Eval took: ", time_eval / 60, " minutes.")

			if not isinstance(metrics, pd.DataFrame):
				metrics = pd.DataFrame([metrics])

			metrics_path = "results/metrics/"

			metrics.to_csv(
				metrics_path + f"{basename}_metrics.csv",
				index=False,
			)

			# can add some italy_aq eval here later
			# --------------------------------------

			del model, training_results

	# after you've used up the datset with a set number of replicates, you can delete it
	del data_dict, train_df, val_df, test_df, data_config
	del train_loader, val_loader, test_loader
	torch.cuda.empty_cache()
	import gc

	gc.collect()
