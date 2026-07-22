PLINK2_LINUX_URL = https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250627.zip
PLINK2_DARWIN_URL = https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_arm64_20250627.zip
PLINK1_LINUX_URL = https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20240818.zip
PLINK1_DARWIN_URL = https://s3.amazonaws.com/plink1-assets/plink_mac_20201019.zip

PLINK2 = bin/plink2
PLINK1 = bin/plink

UNAME_S := $(shell uname -s)
ARCH := $(shell uname -m)

.PHONY: all clean install_plink1 install_plink2 install 
.EXPORT_ALL_VARIABLES:

all:
	$(MAKE) install
	$(MAKE) install_plink1
	$(MAKE) install_plink2

activate:
	source ./.venv/bin/activate

install-dev:
	@echo "Installing development dependencies..."
	uv sync --extra dev

install-test:
	@echo "Installing test dependencies..."
	uv sync --extra test

install-all:
	@echo "Installing all dependencies..."
	uv sync --extra test --extra dev
	uv run pre-commit install

install:
	$(MAKE) install-all

install_plink1:
	@if [ "$(UNAME_S)" = "Linux" ]; then \
		echo "Detected Linux $(ARCH)"; \
		echo "Downloading and installing PLINK 1.9..."; \
		curl -LO $(PLINK1_LINUX_URL); \
		unzip plink_linux_x86_64_20240818.zip; \
		mv plink $(PLINK1); \
		chmod +x $(PLINK1); \
		rm plink_linux_x86_64_20240818.zip; \
		echo "PLINK 1.9 installation completed."; \
	elif [ "$(UNAME_S)" = "Darwin" ]; then \
		echo "Detected macOS $(ARCH)"; \
		echo "Downloading and installing PLINK 1.9..."; \
		curl -LO $(PLINK1_DARWIN_URL); \
		unzip plink_mac_20201019.zip; \
		mv plink $(PLINK1); \
		chmod +x $(PLINK1); \
		rm plink_mac_20201019.zip; \
		echo "PLINK 1.9 installation completed."; \
	else \
		echo "Unsupported platform: $(UNAME_S) $(ARCH)"; \
		exit 1; \
	fi

install_plink2:
	@if [ "$(UNAME_S)" = "Linux" ]; then \
		echo "Detected Linux $(ARCH)"; \
		echo "Downloading and installing PLINK 2.0..."; \
		curl -LO $(PLINK2_LINUX_URL); \
		unzip plink2_linux_x86_64_20250627.zip; \
		mv plink2 $(PLINK2); \
		chmod +x $(PLINK2); \
		rm plink2_linux_x86_64_20250627.zip; \
		rm -f prettify; \
		rm -f toy.map; \
		rm -f toy.ped; \
		rm -f LICENSE; \
		echo "PLINK 2.0 installation completed."; \
	elif [ "$(UNAME_S)" = "Darwin" ]; then \
		echo "Detected macOS $(ARCH)"; \
		echo "Downloading and installing PLINK 2.0..."; \
		curl -LO $(PLINK2_DARWIN_URL) || { echo "Download failed"; exit 1; }; \
		unzip plink2_mac_arm64_20250627.zip || { echo "Unzip failed"; exit 1; }; \
		mv plink2 $(PLINK2); \
		chmod +x $(PLINK2); \
		rm plink2_mac_arm64_20250627.zip; \
		rm -f prettify; \
		rm -f toy.map; \
		rm -f toy.ped; \
		rm -f LICENSE; \
		echo "PLINK 2.0 installation completed."; \
	else \
		echo "Unsupported platform: $(UNAME_S) $(ARCH)"; \
		exit 1; \
	fi

clean:
	@echo "Cleaning up..."
	@rm -f $(PLINK2)
	@echo "Clean up completed."

build-image:
	docker pull ubuntu