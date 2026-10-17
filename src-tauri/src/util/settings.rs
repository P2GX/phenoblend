//! Module to persist settings including the location of the hp.json file
//!

use dirs::home_dir;
use serde::{Serialize, Deserialize};
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};


/// Settings to persist between sessions.

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PhenoblendSettings {
    hp_json_file: Option<String>,
    hpoa_file: Option<String>,
    gene_to_disease_file: Option<String>
}

impl PhenoblendSettings {
    
    fn empty() -> Self {
        Self {
            hp_json_file: None,
            hpoa_file: None,
            gene_to_disease_file: None,
        }
    }

    pub fn set_hp_json_path(&mut self, hp_json: &str) -> Result<(), String> {
        let path = Path::new(hp_json);
        if ! path.is_file() {
            return Err(format!("Did not find file at {hp_json}"));
        }
        self.hp_json_file = Some(hp_json.to_string());
        println!("set_hp_json_path to {}", hp_json);
        self.save_settings()?;
        Ok(())
    }

    pub fn get_hp_json_path(&self) -> Result<String, String> {
        match &self.hp_json_file {
            Some (hp_json) => Ok(hp_json.clone()),
            None => Err(format!("hp.json file not initialized"))
        }
    }


    pub fn set_hpoa_path(&mut self, hpoa_path: &str) -> Result<(), String> {
        let path = Path::new(hpoa_path);
        if ! path.is_file() {
            return Err(format!("Did not find HPOA file at {hpoa_path}"));
        }
        self.hpoa_file = Some(hpoa_path.to_string());

        self.save_settings()?;
        Ok(())
    }

    pub fn get_hpoa_path(&self) -> Result<String, String> {
        match &self.hpoa_file {
            Some (hpoa) => Ok(hpoa.clone()),
            None => Err(format!("phenotype.hpoa file not initialized"))
        }
    }
     pub fn set_g2d_path(&mut self, g2d_path: &str) -> Result<(), String> {
        let path = Path::new(g2d_path);
        if ! path.is_file() {
            return Err(format!("Did not find gene-to-disease file at {g2d_path}"));
        }
        self.gene_to_disease_file = Some(g2d_path.to_string());

        
        self.save_settings()?;
        Ok(())
    }

    pub fn get_g2d_path(&self) -> Result<String, String> {
        match &self.gene_to_disease_file {
            Some (g2d) => Ok(g2d.clone()),
            None => Err(format!("gene_to_disease file not initialized"))
        }
    }


   


    pub fn load_settings() -> PhenoblendSettings {
        let _ = ensure_config_directory();
        let path = get_config_file();
        if path.is_err() {
            return PhenoblendSettings::empty();
        }
        let path = path.unwrap();
    
        if !path.exists() {
            // Write default settings if file doesn't exist
            let default_settings = PhenoblendSettings::empty();
            if let Ok(toml_string) = toml::to_string_pretty(&default_settings) {
                let _ = fs::write(&path, toml_string);
            }
            return default_settings;
        }

        fs::read_to_string(&path)
            .and_then(|contents| toml::from_str(&contents).map_err(|e| {
                std::io::Error::new(std::io::ErrorKind::InvalidData, e)
            }))
            .unwrap_or_else(|err| {
                eprintln!("Warning: Failed to load settings.toml: {}. Using defaults.", err);
                PhenoblendSettings::empty()
            })
    }

    pub fn save_settings(&self) -> Result<(), String> {
        let config_file = get_config_file()?;
        let toml_string = toml::to_string_pretty(&self)
            .map_err(|e| format!("Could not serialize settings: {}", e))?;
        println!("Sainv settings toml = {:?}", toml_string);
        let mut file = File::create(config_file)
            .map_err(|e| format!("Could not create settings file: {}", e))?;

        file.write_all(toml_string.as_bytes())
            .map_err(|e| format!("Could not write to settings file: {}", e))?;
        Ok(())

    }


}



fn get_config_path() -> Result<PathBuf, String> {
    match home_dir() {
        Some(mut home) => {
            home.push(".phenoblend");
            Ok(home)
        }
        None => Err(format!("Could not determine Phenoblend home directory"))
    }
}

fn get_config_file() -> Result<PathBuf, String> {
    let mut config_file = get_config_path()?;
    config_file.push("settings.toml"); // ~/.phenoblend/settings.toml
    Ok(config_file)
}

fn ensure_config_directory() -> Result<(), String> {
    let config_dir = get_config_path()?;
    if !config_dir.exists() {
        fs::create_dir_all(&config_dir).expect("Failed to create config directory");
    }
    Ok(())
}


