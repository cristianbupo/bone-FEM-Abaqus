import sys
sys.path.append('../bone-FEM-Abaqus')

import FreeCADPath
import sys
import FreeCAD as App
import geomdl2gmsh as g2g
import FreeCAD2gmsh as f2g
import longBone as lb
import matplotlib.pyplot as plt
import numpy as np


class AngleDimension:
    def __init__(self, center, point1, point2, radius=None, extra_radius=0.2, text="", color="black"):
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

        ax.plot(arc_x, arc_y, color=self.color, linestyle="-", linewidth=0.5)

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
        # angle = np.degrees(mid_angle) - 90
        ax.text(text_x + 0.25, text_y, self.text, color=self.color, ha="center", va="center", fontsize=12) # , rotation=angle)

        radial1_x = cx + final_radius * np.cos(angle1)
        radial1_y = cy + final_radius * np.sin(angle1)
        radial2_x = cx + final_radius * np.cos(angle2)
        radial2_y = cy + final_radius * np.sin(angle2)

        ax.plot([radial1_x, x1], [radial1_y, y1], linestyle="-", color=self.color, linewidth=0.5)
        ax.plot([radial2_x, x2], [radial2_y, y2], linestyle="-", color=self.color, linewidth=0.5)


class Dimension:
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
            ax.plot([x1, offset_x], [y1, y1], linestyle="-", color="black", linewidth=0.5)
            ax.plot([x2, offset_x], [y2, y2], linestyle="-", color="black", linewidth=0.5)

            # Draw the dimension line with arrows
            ax.annotate(
                "",
                xy=(offset_x, y1),
                xytext=(offset_x, y2),
                arrowprops=dict(arrowstyle="<->", color=self.color, linewidth=0.5)
            )

            # Add dimension text
            ax.text(offset_x - 0.4, mid_y, self.text, color=self.color, ha="center", va="bottom", rotation=90, fontsize=12)

        elif self.dimension == "horizontal":
            mid_x = (x1 + x2) / 2
            offset_y = max(y1, y2) + self.offset

            # Draw the extension lines
            ax.plot([x1, x1], [y1, offset_y], linestyle="-", color="black", linewidth=0.5)
            ax.plot([x2, x2], [y2, offset_y], linestyle="-", color="black", linewidth=0.5)

            # Draw the dimension line with arrows
            ax.annotate(
                "",
                xy=(x1, offset_y),
                xytext=(x2, offset_y),
                arrowprops=dict(arrowstyle="<->", color=self.color, linewidth=0.5)
            )

            # Add dimension text
            ax.text(mid_x, offset_y + 0.1, self.text, color=self.color, ha="center", va="bottom", fontsize=12)

        else:
            raise ValueError("Invalid dimension type. Use 'horizontal', 'vertical', or 'auto'.")


fileName='CADs/longBoneDraw.FCStd'
doc = App.openDocument(fileName)
sketch = doc.getObject('Sketch')

sketch = App.ActiveDocument.getObject('Sketch')

bone, _, boneConfig = lb.getBoneData('loadCases/longbone.json')

curves0 = f2g.getBSplineGeom(sketch, 1000)

# fig, ax = g2g.drawContainer(curves0, controlPol=True)

curvesDraw = g2g.extractCurves(curves0, [0, 1, 7])

fig, ax = g2g.drawContainer(curvesDraw, controlPol=True, knotEval=True)

plt.show()

bigCurve = curvesDraw[0]
bigCurveCtrl = bigCurve.ctrlpts
bigCurvePoints = [(pt[0], pt[1]) for pt in bigCurveCtrl]
bigCurveMid = (bigCurve.evaluate_single(0.5)[0], bigCurve.evaluate_single(0.5)[1])

smallCurve = curvesDraw[1]
smallCurveCtrl = smallCurve.ctrlpts
smallCurvePoints = [(pt[0], pt[1]) for pt in smallCurveCtrl]
smallCurveMid = (smallCurve.evaluate_single(0.5)[0], smallCurve.evaluate_single(0.5)[1])

dimension = Dimension(bigCurvePoints[0], bigCurvePoints[-1], offset=-1.0, text=r"$w_b$")
dimension.draw(ax)

dimension = Dimension(smallCurveMid, bigCurvePoints[0], offset=-3.5, text=r"$h_h$")
dimension.draw(ax)

dimension = Dimension(bigCurveMid, smallCurveMid, offset=-3.5, text=r"$t_c$")
dimension.draw(ax)

dimension = Dimension(bigCurvePoints[5], bigCurvePoints[7], offset=-0.5, text=r"$r_x$")
dimension.draw(ax)

dimension = Dimension(bigCurvePoints[7], bigCurvePoints[9], offset=0.5, text=r"$r_y$")
dimension.draw(ax)

center = bigCurvePoints[7]
point1 =  (bigCurvePoints[7][0]+0.1, bigCurvePoints[7][1])
point2 = bigCurvePoints[9]

angleDimension = AngleDimension(center, point1, point2, radius=3.5, text=r"$\alpha_h$")
angleDimension.draw(ax)

center = smallCurvePoints[0]
point1 =  (smallCurvePoints[0][0]+0.1, smallCurvePoints[0][1])
point2 = smallCurvePoints[1]

angleDimension = AngleDimension(center, point1, point2, radius=3.5, text=r"$\alpha_f$")
angleDimension.draw(ax)

# Mark the origin
ax.plot(0, 0, 'ro')  # Red dot at the origin

# Draw the x-axis arrow
ax.arrow(0, 0, 0.5, 0, head_width=0.1, head_length=0.1, fc='black', ec='black')
ax.annotate(r'$x$', xy=(0.6, 0), xytext=(5, -10), textcoords='offset points', fontsize=12, color='black')

# Draw the y-axis arrow
ax.arrow(0, 0, 0, 0.5, head_width=0.1, head_length=0.1, fc='black', ec='black')
ax.annotate(r'$y$', xy=(0, 0.6), xytext=(-10, 5), textcoords='offset points', fontsize=12, color='black')

# Save to svg
# fig.savefig('/home/cristian/git/mainThesisDocument/images/materialsAndMethods/boneCADmodel_noPol.svg',
#             format='svg', dpi=1200) #, transparent=True) Delete background on inkscape

fig, ax = g2g.drawContainer(curves0, controlPol=False)

# Save to svg
# fig.savefig('/home/cristian/git/mainThesisDocument/images/materialsAndMethods/boneCADmodel_noPol.svg', format='svg', dpi=1200) #, transparent=True) Delete background on inkscape
