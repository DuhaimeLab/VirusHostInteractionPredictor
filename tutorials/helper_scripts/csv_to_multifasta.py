'''Helper file to process output of CRISPRCasFinder.'''

import csv


# This function will read the csv file and create an output multifasta file
def csv_to_multifasta(csv_file: str):
    '''Convert csv file into multifasta file.'''
    # Create output filename
    output_filename = ""
    for c in csv_file:
        if c == ".":
            output_filename += ".fasta"
            break
        output_filename += c

    # Create output file stream
    output = open(output_filename, "w", 1)
    print(output_filename)

    # Open the csv file
    with open(csv_file) as csvfile:
        # Careful not to use the filenme as first param in reader
        filereader = csv.reader(csvfile, delimiter=",")

        # Dummy to avoid sending first row (description) to output
        entry = 0

        # Iterate through each row in the csv file
        for row in filereader:
            # Update entry number
            entry += 1

            # Sample name is row[0]
            sample = row[0]

            # Sequence id is row[1]
            seq_id = row[1]

            # Spacers/DR filename is row[2], but not needed.
            # File contents is row[4]
            content = row[3]

            # Create a string to work with
            info = ""

            # Space count for multiple sequences in content
            spaces = 0

            # Iterate character by character in content
            for c in content:
                # Add char to info
                info += c

                # Check if it is the beginning of a header
                if c == ">":
                    # Add sample and id
                    info += sample
                    info += " "
                    info += seq_id
                    info += " "

                # Check if character is a blank space, signifying a new line
                if c == " ":
                    spaces += 1
                    if spaces > 1:
                        info += "\n"
                        spaces = 0
                    info += "\n"

            # Now send info to output file
            if entry > 1:
                output.write(info)
                output.write("\n\n")

        # end for
    # end with


###################### MAIN ####################################

string = input("Enter csv file: ")
print("Concatenating to multifasta...")
csv_to_multifasta(string)
print("Done! Output has .fasta extension")
