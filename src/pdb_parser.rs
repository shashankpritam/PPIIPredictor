use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Structure {
    models: Vec<Model>,
}

impl Structure {
    pub fn new() -> Self {
        Structure {
            models: Vec::new(),
        }
    }
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.serial_number == other.serial_number
            && self.chains == other.chains
    }
}

impl Eq for Model {}

#[derive(Debug, Clone)]
pub struct Model {
    pub serial_number: isize,
    pub chains: Vec<Chain>,
}

#[derive(Debug, Clone)]
pub struct Chain {
    pub id: char,
    pub residues: Vec<Residue>,
}

#[derive(Debug, Clone)]
pub struct Residue {
    pub name: String,
    pub id: isize,
    pub atoms: Vec<Atom>,
}

#[derive(Debug, Clone, Eq, PartialEq, PartialOrd)]
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

pub fn parse_pdb_file(file_path: &str) -> Result<Structure, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
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

    for line in reader.lines() {
        let record_type = &line?.chars().take(6).collect::<String>();

        match record_type.as_ref() {
            "MODEL " => {
                current_model = Model {
                    serial_number: line?
                        .chars()
                        .skip(10)
                        .take(4)
                        .collect::<String>()
                        .trim()
                        .parse::<isize>()?,
                    chains: Vec::new(),
                };
            }
            "ATOM  " | "HETATM" => {
    let atom = Atom {
        serial: line?.chars().skip(6).take(5).collect::<String>().trim().parse()?,
        name: line?.chars().skip(12).take(4).collect(),
        alt_loc: line?.chars().nth(16).unwrap_or(' '),
        res_name: line?.chars().skip(17).take(3).collect(),
        chain_id: line?.chars().nth(21).unwrap_or(' '),
        res_seq: line?.chars().skip(22).take(4).collect::<String>().trim().parse()?,
        icode: line?.chars().nth(26).unwrap_or(' '),
        x: line?.chars().skip(30).take(8).collect::<String>().trim().parse()?,
        y: line?.chars().skip(38).take(8).collect::<String>().trim().parse()?,
        z: line?.chars().skip(46).take(8).collect::<String>().trim().parse()?,
        occupancy: line?.chars().skip(54).take(6).collect::<String>().trim().parse()?,
        temp_factor: line?.chars().skip(60).take(6).collect::<String>().trim().parse()?,
        element: line?.chars().skip(76).take(2).collect(),
        charge: line?.chars().skip(78).take(2).collect(),
    };
    if let Some(index) = current_chain.residues.iter().position(|r| r.id == current_residue.id) {
        // residue exists, add atom to it
        current_chain.residues[index].atoms.push(atom);
    } else {
        // residue does not exist, create new residue and add atom to it
        current_residue.atoms.push(atom);
        current_chain.residues.push(current_residue.clone());
    }
    current_residue.id = atom.res_seq;
    current_residue.name = atom.res_name.clone();
},
"TER   " => {
    current_chain.residues.push(current_residue.clone());
    current_model.chains.push(current_chain.clone());
    current_residue = Residue {
        name: String::new(),
        id: 0,
        atoms: Vec::new(),
    };
    current_chain = Chain {
        id: ' ',
        residues: Vec::new(),
    };
},
"ENDMDL" => {
    current_chain.residues.push(current_residue.clone());
    current_model.chains.push(current_chain.clone());
    structure.models.push(current_model.clone());
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
},
_ => (),
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