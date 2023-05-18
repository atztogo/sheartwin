"""Command-line utility of sheartwin."""
import click


@click.group()
@click.help_option("-h", "--help")
def cmd_root():
    """Command-line utility of sheartwin."""
    pass
