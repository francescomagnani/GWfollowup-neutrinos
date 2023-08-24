import argparse
import tomli

from GUI.tkinterGUI import tkGUI


if __name__ == "__main__":
    print("Please wait, while the application is being launched...")
    print("\n")
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, help="Path to configuration file.")
    parser.add_argument("--iteration", default="0", type=int)

    args = parser.parse_args()
    iteratio = args.iteration


    with open(args.config, mode="rb") as fc:
        configu = tomli.load(fc)

    # graphic interface
    tkGUI(configu,iteratio)
