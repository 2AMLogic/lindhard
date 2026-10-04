use clap::Parser;

/// Monte Carlo transport of ions and electrons in matter.
#[derive(Parser)]
#[command(name = "lindhard", version = lindhard::VERSION, about)]
struct Cli {}

fn main() -> anyhow::Result<()> {
    let _cli = Cli::parse();
    eprintln!(
        "lindhard {}: no simulation engine yet — see docs/architecture.md",
        lindhard::VERSION
    );
    Ok(())
}
