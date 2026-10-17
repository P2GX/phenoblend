# Installation


Phenoblend is available as prepackaged installers for macOS, Windows, and Linux. Download the latest version from the [Releases](https://github.com/p2gx/phenoblend/releases) page.

## Installing on macOS

> **File to download:** `phenoblend_<version>_aarch64.dmg`  
> This is the macOS installer for Apple Silicon (M1/M2/M3/M4 Macs)

Because this application is open-source and distributed for free, it is not signed or notarized by Apple. macOS will warn you the first time you try to open it. Here's how to install:

1. Download the `.dmg` file from the [Releases](https://github.com/p2gx/phenoblend/releases) page
2. Open the DMG and drag the app into your Applications folder
3. When you try to open it, macOS will show the message:  
   *"'phenoblend' is an app downloaded from the Internet. Are you sure you want to open it?"*
 
If you confirm, the app will be started. From now one you can open phenoblend like any other application (e.g., using the Spotlight searchbar).


## Installing on Windows

> **File to download:** `phenoblend_<version>_x64_en-US.msi`  
> Windows installer (MSI format)

1. Download the `.msi` installer from the [Releases](https://github.com/p2gx/phenoblend/releases) page
2. Double-click to start the installer
3. If Windows shows a blue SmartScreen dialog saying:  
   *"Windows protected your PC"*
4. Click **"More info"** → **"Run anyway"**

> **Note:** Windows shows this for unsigned apps from new developers. Once you install and run it, the warning will not reappear.

## Installing on Linux

### Debian/Ubuntu (recommended)

> **File to download:** `phenoblend_<version>_amd64.deb`  
> Debian/Ubuntu package

1. Download the `.deb` package from the [Releases](https://github.com/p2gx/phenoblend/releases) page
2. Install using:
```bash
sudo apt install ./phenoblend_<version>_amd64.deb
```

Or using dpkg:
```bash
sudo dpkg -i phenoblend_<version>_amd64.deb
```

### Other Linux Distributions

> **File to download:** `phenoblend_<version>_amd64.AppImage`  
> Universal Linux application (no installation needed)

1. Download the `.AppImage` file from the [Releases](https://github.com/p2gx/phenoblend/releases) page
2. Make it executable:
```bash
chmod +x phenoblend_<version>_amd64.AppImage
```

3. Run it:
```bash
./phenoblend_<version>_amd64.AppImage
```

**Note for Arch-based distributions (e.g. EndeavourOS, Manjaro, etc.):** The AppImage may not build/run correctly due to `linuxdeploy` compatibility issues. Use the `.deb` package instead (installable via `debtap`), or build from source and run the binary directly at `./src-tauri/target/release/phenoblend`.

### Troubleshooting on Linux

#### Crash on launch: `Failed to create GBM buffer` (hybrid GPU / Wayland)

It seems that the Tauri/wry rendering stack on Linux goes through WebKit2GTK, which by default attempts GPU-accelerated compositing via EGL and GBM. On hybrid GPU systems with NVIDIA drivers under Wayland, GBM buffer allocation can fail because the NVIDIA proprietary driver's GBM support is incomplete or incompatible with how WebKit2GTK requests buffers. The workaround is to disable WebKit2GTK GPU compositing:

```bash
WEBKIT_DISABLE_COMPOSITING_MODE=1 ./src-tauri/target/release/phenoblend
```
