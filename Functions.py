# Funkcia precita jellyfish histo file a vrati dictionary histo

def read_histo_from_file(filename):
    histo = {}
    with open(filename, newline='\n') as f:
        for line in f:
            parsed_line = [int(s) for s in line.strip().split()]
            if len(parsed_line) != 2:
                print("WARNING: parsed line from histo txt looks like this: " + str(parsed_line))
                continue
            nasobnost, pocet = parsed_line[:2]
            histo[nasobnost] = pocet
    return histo
