import matplotlib.pyplot as plt

class Dimension:
    def __init__(self, point1, point2, offset=0.5, text="", color="blue"):
        """
        Initializes a Dimension object.

        Parameters:
            point1: Tuple (x1, y1) - Start point of the measured line.
            point2: Tuple (x2, y2) - End point of the measured line.
            offset: float - Distance of the dimension line from the measured line.
            text: str - Text to display on the dimension line.
            color: str - Color of the dimension line and text.
        """
        self.point1 = point1
        self.point2 = point2
        self.offset = offset
        self.text = text
        self.color = color

    def draw(self, ax):
        """
        Draws the dimension on the given Matplotlib axis.

        Parameters:
            ax: Matplotlib axis - The axis to draw the dimension on.
        """
        x1, y1 = self.point1
        x2, y2 = self.point2
        dx, dy = x2 - x1, y2 - y1

        # Calculate perpendicular direction for offset
        perp_dx, perp_dy = -dy, dx
        perp_length = (perp_dx**2 + perp_dy**2)**0.5
        perp_dx, perp_dy = perp_dx / perp_length, perp_dy / perp_length

        # Offset points for the dimension line
        offset_start = (x1 + self.offset * perp_dx, y1 + self.offset * perp_dy)
        offset_end = (x2 + self.offset * perp_dx, y2 + self.offset * perp_dy)

        # Draw the extension lines
        ax.plot([x1, offset_start[0]], [y1, offset_start[1]], linestyle="--", color="black")
        ax.plot([x2, offset_end[0]], [y2, offset_end[1]], linestyle="--", color="black")

        # Draw the dimension line with arrows
        ax.annotate(
            "",
            xy=offset_start,
            xytext=offset_end,
            arrowprops=dict(arrowstyle="<->", color=self.color)
        )

        # Add dimension text
        mid_x = (offset_start[0] + offset_end[0]) / 2
        mid_y = (offset_start[1] + offset_end[1]) / 2
        ax.text(mid_x, mid_y, self.text, color=self.color, ha="center", va="bottom")

# Example Usage
fig, ax = plt.subplots()

# Create two points
point1 = (1, 1)
point2 = (5, 3)

# Create a Dimension object
dimension = Dimension(point1, point2, offset=0.5, text="4.47 units", color="red")

# Draw the measured line
ax.plot([point1[0], point2[0]], [point1[1], point2[1]], color="black", linewidth=2)

# Draw the dimension
dimension.draw(ax)

# Adjust the plot
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_aspect("equal")
plt.grid(True)
plt.show()
