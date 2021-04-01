from numpy import append


class ColorSorter:
    def __init__(self, f_colors, node_order, no_colors):
        self.f_colors = f_colors
        self.node_order = node_order
        self.N = no_colors
        self.pos = {node_order[i]: i for i in range(len(node_order))}

    def rw(self, c_A, c_B):
        """Calculates a 'relative weight' of the color."""
        for a in c_A:
            for b in c_B:
                if self.pos[a] > self.pos[b]:
                    RW += 1
                elif self.pos[a] < self.pos[b]:
                    RW -= 1
        return RW

    def sort_colors(self):
        colors = [c for c in range(self.N)]
        changed = True

        while changed:
            changed = False
            for i in range(self.N-1:
                if self.rw(colors[i], colors[i+1]) > 0:
                    c = colors[i]
                    colors[i] = colors[i+1]
                    colors[i] = c
                    changed = True

        customers = {c: [] for c in range(no_colors)}

        for j in node_order:
            customers[f_colors[j]].append(j)

        return colors, customers
