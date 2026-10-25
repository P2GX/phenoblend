# Developers <Badge type="tip" text="v1.0" />

Welcome to the development guide for **Phenoblend**. Follow the steps below to get your local environment running.

## Getting Started

<Steps>

### 1. Clone the Repository
Clone the source code from the GitHub repository:

```bash
git clone git@github.com:P2GX/pyphenoblend.git
```

### 2. Install Node.js
If you haven't already, install Node.js and npm by following [the official npm instructions](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm).


### 3. Install Dependencies
Navigate into the project directory and install the required packages:
Change into the directory and install the angular packages

```bash
cd phenoblend
npm install
```

### 4. Run in Development Mode

Start the Tauri development server:


```bash
npm run tauri dev
```

::: tip Tip for Tauri
Make sure you have all the required system dependencies installed for Tauri (such as Rust and webview development headers) before running the dev command. Check the Tauri Prerequisites Guide if you run into issues.
:::

</Steps>

### Contributing & Support

If you encounter a bug or have a feature request, please let us know by opening an issue on GitHub.

::: warning Found a bug?
Please provide steps to reproduce the issue along with your operating system and environment details.
:::

* [GitHub Issues Tracker](https://github.com/P2GX/phenoblend/issues)



## Building from Source

This will work on any OS.

```bash
git clone https://github.com/P2GX/phenoblend.git
cd phenoblend
npm install
npm run tauri build
```
The built installers will appear under:
```bash
src-tauri/target/release/bundle/
``` 


## Prerequisites

### Node.js (at least version 18) and npm (at least version 9). 

You can check if you have them installed via

```bash
node -v
npm -v
``` 

If necessary, go to [Node.js](https://nodejs.org) to install these programs.

### Rust and Cargo

See [rustup](https://rustup.rs) if needed.

### Git 

If you do not have git installed, replace the cloning step with a download of the archive.

## Platform-specific code
Please report any dependencies not listed above.

# Running the app

1. Clone from GitHub
```bash 
git clone https://github.com/P2GX/phenoblend.git
cd phenoblend
```

2. Install npm dependencies
From within the phenoblend directory, enter
```bash
npm install
```

3. Running the app
```bash
npm run tauri dev
```

This will run the application.


### Documentation

Phenoblend's documentation is created using vitepress. It can be run locally as follows.

```bash
npm run docs:dev
```

Look for a link such as `http://localhost:5173/phenoblend/` in the shell and click on it. Corrections and suggestions for the documentation are also welcome as issues or pull requests.


### Running the batch analysis


```bash
cd src-tauri
cargo run --release --features cli --bin phenoblend-batch -- --phenopackets ../ppkts --hpo hp.json --hpoa phenotype.hpoa --out results.tsv --shuffle-seed 42
```


TODO: details