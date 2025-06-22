import intermol

print(f"InterMol version: {intermol.__version__}")
print("-" * 30)

# See what's directly in the 'intermol' module
print("Available in 'intermol':")
for item in dir(intermol):
    if not item.startswith('_'): # Exclude private/special attributes
        print(f"  {item}")
print("-" * 30)

# Check if 'amber' is a submodule and what it contains
if hasattr(intermol, 'amber'):
    print("Available in 'intermol.amber':")
    for item in dir(intermol.amber):
        if not item.startswith('_'):
            print(f"  intermol.amber.{item}")
    print("-" * 30)
    # If AmberParser was in the previous attempt, maybe it's directly under intermol.amber
    if hasattr(intermol.amber, 'AmberParser'):
        print("Found intermol.amber.AmberParser!")
    else:
        print("intermol.amber.AmberParser NOT found directly.")
else:
    print("'intermol.amber' submodule NOT found.")
print("-" * 30)

# Check if 'gromacs' is a submodule and what it contains
if hasattr(intermol, 'gromacs'):
    print("Available in 'intermol.gromacs':")
    for item in dir(intermol.gromacs):
        if not item.startswith('_'):
            print(f"  intermol.gromacs.{item}")
    print("-" * 30)
    if hasattr(intermol.gromacs, 'GromacsDriver'):
         print("Found intermol.gromacs.GromacsDriver!")
    elif hasattr(intermol.gromacs, 'GromacsWriter'): # Common alternative name
         print("Found intermol.gromacs.GromacsWriter!")
    else:
        print("intermol.gromacs.GromacsDriver/Writer NOT found directly.")
else:
    print("'intermol.gromacs' submodule NOT found.")
print("-" * 30)

# The rest of the script is commented out as we first need to find the correct API calls
"""
from intermol.system import System 

print("Starting InterMol conversion...")

# Define input files
prmtop_file = 'mixture.prmtop'
coord_file = 'mixture.inpcrd'

# Define output files
gromacs_top_out = 'mixture_intermol.top'
gromacs_gro_out = 'mixture_intermol.gro'

current_system = System(name="mixture_system")

try:
    # Based on the output of dir() above, we would fill in the correct
    # way to initialize the Amber parser and read files, e.g.:
    #
    # if hasattr(intermol.amber, 'AmberParser'):
    #     amber_parser = intermol.amber.AmberParser(prmtop_filename=prmtop_file, crd_filename=coord_file, system_name="mixture_system")
    #     amber_parser.read()
    #     system_from_amber = amber_parser.system 
    #     print("Amber files loaded using AmberParser.")
    # else:
    #     raise NotImplementedError("Could not find a way to parse Amber files with this InterMol version via script.")

    # And for writing:
    # if hasattr(intermol.gromacs, 'GromacsDriver'):
    #     gromacs_driver = intermol.gromacs.GromacsDriver()
    #     gromacs_driver.write(system_from_amber, top=gromacs_top_out, gro=gromacs_gro_out)
    #     print("GROMACS files written using GromacsDriver.")
    # elif hasattr(intermol.gromacs, 'GromacsWriter'):
    #      gromacs_writer = intermol.gromacs.GromacsWriter() # Or however it's initialized
    #      gromacs_writer.write_files(system_from_amber, top=gromacs_top_out, gro=gromacs_gro_out) # Method name might vary
    #      print("GROMACS files written using GromacsWriter.")
    # else:
    #      raise NotImplementedError("Could not find a way to write GROMACS files with this InterMol version via script.")

    print("Conversion using InterMol done.")

except AttributeError as e:
    print(f"AttributeError during InterMol operations: {e}")
except ImportError as e:
    print(f"ImportError during InterMol operations: {e}")
except Exception as e:
    print(f"An unexpected error occurred during InterMol operations: {e}")
"""
print("Script finished exploration. Next step is to use the printed API info.")
