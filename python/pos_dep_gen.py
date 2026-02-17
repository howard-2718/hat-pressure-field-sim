# Generates an integer lattice spiral, centered at (0, 0). 
# The starting direction can be chosed, as well as the stopping condition.
# When the points generated from this spiral are converted into positions in the pressure field, 
def generate_spiral():
    # Up, Right, Down, Left
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]

    x = y = 0

    step_len = 1