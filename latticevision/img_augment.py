import random
import torch


def random_negative_clim(fields, params):
	"""
	Randomly negates fields, but has
	no impact on parameters.
	"""
	fields = -fields
	return fields, params


def random_translate_clim(fields, params):
	"""
	Randomly translates a pair of fields and parameters
	by a random shift in the x and y directions.
	"""
	# get dims for bounding the shift
	_, _, rows, cols = fields.shape

	# randomly shift (at max shift whole image)
	xshift = random.randint(-cols, cols)
	yshift = random.randint(-rows, rows)

	# roll fields
	fields = torch.roll(fields, shifts=xshift, dims=-1)
	fields = torch.roll(fields, shifts=yshift, dims=-2)

	# roll second field the same way (param)
	params = torch.roll(params, shifts=xshift, dims=-1)
	params = torch.roll(params, shifts=yshift, dims=-2)

	return fields, params


def random_augment_clim(fields, params):
	"""
	Randomly augments a pair of fields and parameters
	by applying one or more augmentations.
	"""

	# available augmentations
	augs = [random_negative_clim, random_translate_clim]
	# randomly choose one or more of them
	num_augs = random.randint(1, len(augs))
	chosen_augs = random.sample(augs, num_augs)

	# apply to both fields and params
	for aug in chosen_augs:
		fields, params = aug(fields, params)

	return fields, params
