# gerber_to_gcode
Simple python script for converting gerber files into gcode to directly print solder paste. It depends mainly on the project [gerber_to_scad](https://github.com/kirberich/gerber_to_scad).

## Installation

```bash
cd gerber_to_gcode
virtualenv env # Only tested for python 2.7
source env/bin/activate
pip install -r requirements.txt
```

You should now be able to run the script. You'll get some information on available options if you run it with the -h argument:
```bash
(env) $ python gerber_to_gcode.py -h
usage: gerber_to_gcode.py [-h] [-t THICKNESS] [-x OFFSET_X] [-y OFFSET_Y]
                          [-z OFFSET_Z] [-s SOLDER_HEIGHT] [-f FLOWRATE]
                          [-r RETRACTION]
                          outline_file solderpaste_file output_file

Convert gerber files to gcode to print solderpaste with a 3d printer.

positional arguments:
  outline_file          Outline file
  solderpaste_file      Solderpaste file
  output_file           Output file

optional arguments:
  -h, --help            show this help message and exit
  -t THICKNESS, --thickness THICKNESS
                        Thickness (in mm) of the PCB. (default: 1.6)
  -x OFFSET_X, --offset_x OFFSET_X
                        Offset in X-direction. (default: 0.0)
  -y OFFSET_Y, --offset_y OFFSET_Y
                        Offset in Y-direction. (default: 0.0)
  -z OFFSET_Z, --offset_z OFFSET_Z
                        Offset in Z-direction. (default: 0.0)
  -s SOLDER_HEIGHT, --solder_height SOLDER_HEIGHT
                        Height of the solder paste. (default: 0.3)
  -f FLOWRATE, --flowrate FLOWRATE
                        Increase the flow rate (in %) of the solder paste.
                        (default: 100)
  -r RETRACTION, --retraction RETRACTION
                        Retraction length (in mm) of the solder paste.
                        (default: 2.0)
```

For basic usage, simply run the script with input files for the gerber outline and solderpaste files and specify an output:

```bash
python gerber_to_gcode.py outline_file.gko toppaste_file.gtp output.gcode
```
