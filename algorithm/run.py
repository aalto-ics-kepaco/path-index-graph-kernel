import argparse

parser = argparse.ArgumentParser(description='Run the path index graph kernel.')
parser.add_argument('-m', '--molfolder', type=str, nargs=1, dest='molfolder', required=True,
                   help='source folder with .mol or .sdf files')

args = parser.parse_args()
print args.molfolder

# @TODO: finish script