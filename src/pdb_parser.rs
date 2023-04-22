use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Structure {
    pub models: Vec<Model>,
    pub chains: Vec<Chain>,
}

impl PartialEq for Structure {
    fn eq(&self, other: &Self) -> bool {
        self.models.len() == other.models.len()
            && self.models.iter().zip(other.models.iter()).all(|(a, b)| a == b)
            && self.chains.len() == other.chains.len()
            && self.chains.iter().zip(other.chains.iter()).all(|(a, b)| a == b)
    }
}

impl Structure {
    pub fn new() -> Self {
        Structure {
            models: Vec::new(),
            chains: Vec::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Model {
    pub serial_number: isize,
    pub chains: Vec<Chain>,
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.serial_number == other.serial_number
            && self.chains.iter().zip(other.chains.iter()).all(|(a, b)| a == b)
    }
}

impl Eq for Model {}

#[derive(Debug, Clone)]
pub struct Chain {
    pub id: char,
    pub residues: Vec<Residue>,
}

impl PartialEq for Chain {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
            && self.residues.len() == other.residues.len()
            && self.residues.iter().zip(other.residues.iter()).all(|(a, b)| a == b)
    }
}


#[derive(Debug, Clone)]
pub struct Residue {
    pub name: String,
    pub id: isize,
    pub atoms: Vec<Atom>,
}

impl PartialEq for Residue {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
            && self.id == other.id
            && self.atoms.iter().zip(other.atoms.iter()).all(|(a, b)| a == b)
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub serial: isize,
    pub name: String,
    pub alt_loc: char,
    pub res_name: String,
    pub chain_id: char,
    pub res_seq: isize,
    pub icode: char,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub occupancy: f32,
    pub temp_factor: f32,
    pub element: String,
    pub charge: String,
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.serial == other.serial
            && self.name == other.name
            && self.alt_loc == other.alt_loc
            && self.res_name == other.res_name
            && self.chain_id == other.chain_id
            && self.res_seq == other.res_seq
            && self.icode == other.icode
            && (self.x - other.x).abs() < f32::EPSILON
            && (self.y - other.y).abs() < f32::EPSILON
            && (self.z - other.z).abs() < f32::EPSILON
            && (self.occupancy - other.occupancy).abs() < f32::EPSILON
            && (self.temp_factor - other.temp_factor).abs() < f32::EPSILON
            && self.element == other.element
            && self.charge == other.charge
    }
}


pub fn parse_pdb_file(file_content: &str) -> Result<Structure, Box<dyn Error>> {
    let reader = BufReader::new(file_content.as_bytes());
    let mut structure = Structure::new();
    let mut current_model = Model {
        serial_number: 1,
        chains: Vec::new(),
    };

    let mut current_chain = Chain {
        id: ' ',
        residues: Vec::new(),
    };

    let mut current_residue = Residue {
        name: String::new(),
        id: 0,
        atoms: Vec::new(),
    };

    for result_line in reader.lines() {
        let line = result_line?;
        let record_type = line.chars().take(6).collect::<String>();
        //println!("Line: {}", line);
        match record_type.as_ref() {
            // A new model is encountered
            "MODEL " => {
                current_model = Model {
                    serial_number: line
                        .chars()
                        .skip(10)
                        .take(4)
                        .collect::<String>()
                        .trim()
                        .parse::<isize>()?,
                    chains: Vec::new(),
                };
                //println!("New model encountered: {:?}", current_model);
            }
            // An atom is encountered
            "ATOM  " | "HETATM" => {
                let atom = Atom {
                    serial: line.chars().skip(6).take(5).collect::<String>().trim().parse()?,
                    name: line.chars().skip(12).take(4).collect(),
                    alt_loc: line.chars().nth(16).unwrap_or(' '),
                    res_name: line.chars().skip(17).take(3).collect(),
                    chain_id: line.chars().nth(21).unwrap_or(' '),
                    res_seq: line.chars().skip(22).take(4).collect::<String>().trim().parse()?,
                    icode: line.chars().nth(26).unwrap_or(' '),
                    x: line.chars().skip(30).take(8).collect::<String>().trim().parse()?,
                    y: line.chars().skip(38).take(8).collect::<String>().trim().parse()?,
                    z: line.chars().skip(46).take(8).collect::<String>().trim().parse()?,
                    occupancy: line.chars().skip(54).take(6).collect::<String>().trim().parse()?,
                    temp_factor: line.chars().skip(60).take(6).collect::<String>().trim().parse()?,
                    element: line.chars().skip(76).take(2).collect(),
                    charge: line.chars().skip(78).take(2).collect(),
                };
                //println!("Atom encountered: {:?}", atom);

                // Check if the current residue exists in the current chain
                if let Some(index) = current_chain.residues.iter().position(|r| r.id == current_residue.id) {
                    // If the residue exists, add the atom to it
                    current_chain.residues[index].atoms.push(atom.clone());
                } else {
                    // If the residue does not exist, create a new residue and add the atom to it
                    current_residue.atoms.push(atom.clone());
                    current_chain.residues.push(current_residue.clone());
                }

                // Update the current residue's ID and name
                current_residue.id = atom.res_seq;
                current_residue.name = atom.res_name.clone();
                //println!("Current residue updated: {:?}", current_residue);
            }
            // A chain terminator is encountered
            "TER   " => {
                current_chain.residues.push(current_residue.clone());
                current_model.chains.push(current_chain.clone());
                //println!("TER encountered: current_chain: {:?}", current_chain);
                current_chain = Chain {
                    id: ' ',
                    residues: Vec::new(),
                };
                //println!("Resetting current_chain: {:?}", current_chain);
            }

            "ENDMDL" => {
                current_chain.residues.push(current_residue.clone());
                current_model.chains.push(current_chain.clone());
                //println!("ENDMDL encountered: current_model: {:?}", current_model);

                // Push the current model to the structure
                structure.models.push(current_model.clone());

                // Reset the current_residue, current_chain, and current_model after pushing them to the structure
                current_residue = Residue {
                    name: String::new(),
                    id: 0,
                    atoms: Vec::new(),
                };
                current_chain = Chain {
                    id: ' ',
                    residues: Vec::new(),
                };
                current_model = Model {
                    serial_number: current_model.serial_number + 1,
                    chains: Vec::new(),
                };
                //println!("Resetting current_residue, current_chain, and current_model");
            }
            _ => (),
        }
    }
    // Check if the current residue has any atoms and add it to the current chain
    if !current_residue.atoms.is_empty() {
        current_chain.residues.push(current_residue.clone());
    }

    // Check if the current chain has any residues and add it to the current model
    if !current_chain.residues.is_empty() {
        current_model.chains.push(current_chain.clone());
    }

    // Check if the current model has any chains and add it to the structure
    if !current_model.chains.is_empty() {
        structure.models.push(current_model.clone());
    }

    Ok(structure)
}





pub fn parse_atom(line: &str) -> Atom {
    let serial: isize = line[6..11].trim().parse().unwrap();
    let name = line[12..16].trim().to_string();
    let alt_loc = line[16..17].chars().next().unwrap_or(' ');
    let res_name = line[17..20].trim().to_string();
    let chain_id = line[21..22].chars().next().unwrap_or(' ');
    let res_seq: isize = line[22..26].trim().parse().unwrap();
    let icode = line[26..27].chars().next().unwrap_or(' ');
    let x: f32 = line[30..38].trim().parse().unwrap();
    let y: f32 = line[38..46].trim().parse().unwrap();
    let z: f32 = line[46..54].trim().parse().unwrap();
    let occupancy: f32 = line[54..60].trim().parse().unwrap_or(0.0);
    let temp_factor: f32 = line[60..66].trim().parse().unwrap_or(0.0);
    let element = line[76..78].trim().to_string();
    let charge = line[78..80].trim().to_string();

    Atom {
        serial,
        name,
        alt_loc,
        res_name,
        chain_id,
        res_seq,
        icode,
        x,
        y,
        z,
        occupancy,
        temp_factor,
        element,
        charge,
    }
}