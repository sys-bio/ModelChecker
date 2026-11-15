//! Converts an Antimony file to SBML, for making test and demo models.
//!
//! cargo run -p mc-core --example ant2sbml -- model.ant > model.xml

fn main() {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("usage: ant2sbml FILE.ant");
        std::process::exit(2);
    };
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("{path}: {e}");
        std::process::exit(1);
    });
    let loaded = antimony::parse(&text);
    for d in loaded.diagnostics.iter() {
        eprintln!("{d}");
    }
    if loaded.has_errors() {
        std::process::exit(1);
    }
    print!("{}", antimony::sbml::write_sbml(&loaded.document));
}
