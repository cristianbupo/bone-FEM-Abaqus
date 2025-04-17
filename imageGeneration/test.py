import matplotlib.pyplot as plt
import numpy as np

class Dimension:
    def __init__(self, point1, point2, offset=0.5, text="", color="black"):
        self.point1 = point1
        self.point2 = point2
        self.offset = offset
        self.text = text
        self.color = color

    def draw(self, ax):
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
        ax.plot([x1, offset_start[0]], [y1, offset_start[1]], linestyle="-", color="black")
        ax.plot([x2, offset_end[0]], [y2, offset_end[1]], linestyle="-", color="black")

        # Draw the dimension line with arrows
        ax.annotate(
            "",
            xy=offset_start,
            xytext=offset_end,
            arrowprops=dict(arrowstyle="<->", color=self.color)
        )

        # Add dimension text with rotation
        mid_x = (offset_start[0] + offset_end[0]) / 2
        mid_y = (offset_start[1] + offset_end[1]) / 2
        angle = np.degrees(np.arctan2(dy, dx))
        ax.text(mid_x, mid_y, self.text, color=self.color, ha="center", va="bottom", rotation=angle)

class AngleDimension:
    def __init__(self, center, point1, point2, radius=None, extra_radius=0.2, text="", color="green"):
        self.center = center
        self.point1 = point1
        self.point2 = point2
        self.radius = radius
        self.extra_radius = extra_radius
        self.text = text
        self.color = color

    def draw(self, ax):
        cx, cy = self.center
        x1, y1 = self.point1
        x2, y2 = self.point2

        dist1 = np.sqrt((x1 - cx)**2 + (y1 - cy)**2)
        dist2 = np.sqrt((x2 - cx)**2 + (y2 - cy)**2)
        base_radius = self.radius if self.radius is not None else max(dist1, dist2)
        final_radius = base_radius + self.extra_radius

        angle1 = np.arctan2(y1 - cy, x1 - cx)
        angle2 = np.arctan2(y2 - cy, x2 - cx)

        if angle1 > angle2:
            angle1, angle2 = angle2, angle1
            x1, y1, x2, y2 = x2, y2, x1, y1

        if angle2 - angle1 > np.pi:
            angle1, angle2 = angle2, angle1 + 2 * np.pi

        theta = np.linspace(angle1, angle2, 100)
        arc_x = cx + final_radius * np.cos(theta)
        arc_y = cy + final_radius * np.sin(theta)

        ax.plot(arc_x, arc_y, color=self.color, linestyle="-")

        arc_start_x = cx + final_radius * np.cos(angle1)
        arc_start_y = cy + final_radius * np.sin(angle1)
        arc_end_x = cx + final_radius * np.cos(angle2)
        arc_end_y = cy + final_radius * np.sin(angle2)

        ax.annotate(
            "",
            xy=(arc_start_x, arc_start_y),
            xytext=(arc_start_x + 0.2 * np.cos(angle1 + np.pi / 2),
                    arc_start_y + 0.2 * np.sin(angle1 + np.pi / 2)),
            arrowprops=dict(arrowstyle="->", color=self.color)
        )

        ax.annotate(
            "",
            xy=(arc_end_x, arc_end_y),
            xytext=(arc_end_x + 0.2 * np.cos(angle2 - np.pi / 2),
                    arc_end_y + 0.2 * np.sin(angle2 - np.pi / 2)),
            arrowprops=dict(arrowstyle="->", color=self.color)
        )

        mid_angle = (angle1 + angle2) / 2
        text_x = cx + (final_radius + 0.2) * np.cos(mid_angle)
        text_y = cy + (final_radius + 0.2) * np.sin(mid_angle)
        angle = np.degrees(mid_angle) - 90
        ax.text(text_x, text_y, self.text, color=self.color, ha="center", va="center", rotation=angle)

        radial1_x = cx + final_radius * np.cos(angle1)
        radial1_y = cy + final_radius * np.sin(angle1)
        radial2_x = cx + final_radius * np.cos(angle2)
        radial2_y = cy + final_radius * np.sin(angle2)

        ax.plot([radial1_x, x1], [radial1_y, y1], linestyle="-", color=self.color)
        ax.plot([radial2_x, x2], [radial2_y, y2], linestyle="-", color=self.color)


class HorizontalVerticalDimension:
    def __init__(self, point1, point2, offset=0.2, text="", color="black", dimension="auto"):
        self.point1 = point1
        self.point2 = point2
        self.offset = offset
        self.text = text
        self.color = color
        self.dimension = dimension

    def draw(self, ax):
        x1, y1 = self.point1
        x2, y2 = self.point2

        if self.dimension == "auto":
            if abs(x2 - x1) > abs(y2 - y1):
                self.dimension = "horizontal"
            else:
                self.dimension = "vertical"

        if self.dimension == "vertical":
            mid_y = (y1 + y2) / 2
            offset_x = max(x1, x2) + self.offset

            # Draw the extension lines
            ax.plot([x1, offset_x], [y1, y1], linestyle="-", color="black")
            ax.plot([x2, offset_x], [y2, y2], linestyle="-", color="black")

            # Draw the dimension line with arrows
            ax.annotate(
                "",
                xy=(offset_x, y1),
                xytext=(offset_x, y2),
                arrowprops=dict(arrowstyle="<->", color=self.color)
            )

            # Add dimension text
            ax.text(offset_x - 0.07, mid_y, self.text, color=self.color, ha="center", va="bottom", rotation=90)

        elif self.dimension == "horizontal":
            mid_x = (x1 + x2) / 2
            offset_y = max(y1, y2) + self.offset

            # Draw the extension lines
            ax.plot([x1, x1], [y1, offset_y], linestyle="-", color="black")
            ax.plot([x2, x2], [y2, offset_y], linestyle="-", color="black")

            # Draw the dimension line with arrows
            ax.annotate(
                "",
                xy=(x1, offset_y),
                xytext=(x2, offset_y),
                arrowprops=dict(arrowstyle="<->", color=self.color)
            )

            # Add dimension text
            ax.text(mid_x, offset_y, self.text, color=self.color, ha="center", va="bottom")

        else:
            raise ValueError("Invalid dimension type. Use 'horizontal', 'vertical', or 'auto'.")


fig, ax = plt.subplots()

ax.set_aspect('equal', adjustable='box')
ax.axis('off')

# Create two points
point1 = (1, 1)
point2 = (5, 3)

# Create a Dimension object
dimension = Dimension(point1, point2, offset=0.5, text=r"$a$", color="black")

# Draw the measured line
# ax.plot([point1[0], point2[0]], [point1[1], point2[1]], color="black", linewidth=2)

# Draw the dimension
# dimension.draw(ax)

# Create two points and a center
center = (3, 3)
point1 = (1, 3)
point2 = (4, 4)

# Create an AngleDimension object
angle_dimension = AngleDimension(center, point1, point2, extra_radius=0.5, text=r"$\alpha$", color="black")

# Draw the angle dimension
# angle_dimension.draw(ax)

# Create two points for horizontal dimension
point1 = (1, 1)
point2 = (5, 3)

# Create a HorizontalVerticalDimension object for horizontal dimension
horizontal_dimension = HorizontalVerticalDimension(point1, point2, offset=0.5, text=r"$h$", color="black", dimension="horizontal")

# Draw the horizontal dimension
# horizontal_dimension.draw(ax)

# Create two points for vertical dimension
point1 = (3, 2)
point2 = (3, 5)

# Create a HorizontalVerticalDimension object for vertical dimension
vertical_dimension = HorizontalVerticalDimension(point1, point2, offset=0.5, text=r"$v$", color="black", dimension="vertical")

# Draw the vertical dimension
# vertical_dimension.draw(ax)

# Create two points for auto dimension
point1 = (2, 2)
point2 = (4, 5)

# Create a HorizontalVerticalDimension object for auto dimension
auto_dimension = HorizontalVerticalDimension(point2, point1, offset=0.5, text=r"$d$", color="black", dimension="vertical")

# Draw the auto dimension
auto_dimension.draw(ax)

plt.legend()
plt.show()