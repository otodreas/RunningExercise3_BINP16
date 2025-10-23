"""
Put this program inside a folder in which you want to remove bad characters
from filenames. Run it from the command line, passing 'R' as the sole argument
if you want to update all the files recursively (all at once).

Oliver Todreas
2025-10-23
"""

# Import libraries
import sys
import os

# Import arguments, set recursive variable
if len(sys.argv) == 1:
    recursive = False
elif len(sys.argv) == 2:
    if sys.argv[1].upper() == "R":
        recursive = True
    else:
        sys.exit("To rename files recursively (all at once): pass 'r'.")
else:
    sys.exit("Too many arguments.")

# Set variables allbadchars, files, and newnames
allbadchars = "#%&{} <>*?/$!:@+`|="
files = os.listdir()
files.sort()
files.pop(files.index(sys.argv[0]))
newnames = {}

# Loop through files, store new names in the dictionary newnames
for file in files:
    if not os.path.isdir(os.path.join(os.getcwd(), file)):  # Do not change folder names
        badchars = list(set(allbadchars) & set(file))
        if len(badchars) > 0:
            for c in badchars:
                file_new = file.replace(c, "_")
                if recursive:
                    newnames[file] = file_new
                else:
                    choice = input(
                        f"Rename {file} -> {file_new}? [N/y] (abort: X) "
                    ).upper()
                    if choice == "Y":
                        newnames[file] = file_new
                    elif choice == "X":
                        sys.exit("Non-recursive rename aborted, no files renamed.")
                    else:
                        print(f"Skipped: {file}.")

# Rename files
for key, value in newnames.items():
    os.rename(key, value)
    print(f"RENAMED: {key} -> {value}.")
