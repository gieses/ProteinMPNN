# change your shell here, this is needed for smooth conda workflows.
SHELL := $(shell which zsh)
.ONESHELL:
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
.DEFAULT_GOAL := help
.PHONY: help conda clone test yml env rm_env rm_yml
ENV_NAME = "proteinmpnn"

help:				## show this help
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e "s/\\$$//" | sed -e "s/##//"

conda:				## install needed conda packages for env creation (not proteinmpnn)
	conda install mamba conda-devenv -y

clone:				## clone the proteinmpnn repository
	git clone https://github.com/gieses/ProteinMPNN

yml:				## create environment.yml using conda devenv
	 conda devenv --file environment.devenv.yml --print >environment.yml

env:				## just create the env from the environment.yml
	mamba env update -f environment.yml

rm_env:				## remove conda env
	conda env remove -n $(ENV_NAME)

rm_yml:				## remove environment.yml
	rm -f environment.yml

pip:
	pip install -e .