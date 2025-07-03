import matplotlib.pyplot as plt
import numpy as np

class RectangleDivider:
# the rectangle is in the baselink frame
# the x-axis is along the length
# the y-axis is along the width
    def __init__(self, vertices):
        self.vertices = vertices
        x_coords = vertices[:, 0]
        y_coords = vertices[:, 1]
        self.vertices = vertices
        self.x1, self.y1 = np.min(x_coords), np.min(y_coords)
        self.x2, self.y2 = np.max(x_coords), np.max(y_coords)
        self.rectangles = [(self.x1, self.y1, self.x2, self.y2)]
        self.circles = []
        self.outer_circles = []
        self.radius = 0
        self.sagitta_lateral = 0
        self.sagitta_longitudinal = 0
        self.previous_max_dimension = max(self.x2 - self.x1, self.y2 - self.y1)

    def divide(self, divisions):
        for _ in range(divisions):
            self._divide_longest_side()
        self._compute_min_covering_circles()

    def _divide_longest_side(self):
        new_rectangles = []
        
        for rect in self.rectangles:
            x1, y1, x2, y2 = rect
            width = x2 - x1
            height = y2 - y1

            if width > height:
                mid_x = (x1 + x2) / 2
                new_rectangles.append((x1, y1, mid_x, y2))
                new_rectangles.append((mid_x, y1, x2, y2))
                new_max_dimension = width / 2
            else:
                mid_y = (y1 + y2) / 2
                new_rectangles.append((x1, y1, x2, mid_y))
                new_rectangles.append((x1, mid_y, x2, y2))
                new_max_dimension = height / 2

        self.previous_max_dimension = new_max_dimension
        self.rectangles = new_rectangles

    def _compute_min_covering_circles(self):
        self.circles = []
        self.outer_circles = []
        for rect in self.rectangles:
            x1, y1, x2, y2 = rect
            center_x = (x1 + x2) / 2
            center_y = (y1 + y2) / 2
            radius = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) / 2
            self.radius = radius
            width = x2 - x1
            height = y2 - y1
            self.sagitta_lateral = radius - np.sqrt(radius**2 - (width / 2)**2)
            self.sagitta_longitudinal = radius - np.sqrt(radius**2 - (height / 2)**2)
            self.circles.append((center_x, center_y))
            if (center_x - radius <= self.x1 or center_x + radius >= self.x2 or
                center_y - radius <= self.y1 or center_y + radius >= self.y2):
                self.outer_circles.append((center_x, center_y))

    def get_rectangles(self):
        return self.rectangles

    def get_circles(self):
        return self.circles

    def get_outer_circles(self):
        return self.outer_circles

    def plot_rectangles_and_circles(self):
        fig, ax = plt.subplots()
        for rect in self.rectangles:
            x1, y1, x2, y2 = rect
            width = x2 - x1
            height = y2 - y1
            ax.add_patch(plt.Rectangle((x1, y1), width, height, edgecolor='blue', facecolor='none', lw=2))
        
        for center_x, center_y in self.outer_circles:
            ax.add_patch(plt.Circle((center_x, center_y), self.radius, edgecolor='red', facecolor='none', lw=1, linestyle='--'))
            
        for center_x, center_y in self.circles:
            if (center_x, center_y) not in self.outer_circles:
                ax.add_patch(plt.Circle((center_x, center_y), self.radius, edgecolor='green', facecolor='none', lw=1))
            
        ax.set_xlim(self.x1 - 1, self.x2 + 1)
        ax.set_ylim(self.y1 - 1, self.y2 + 1)
        ax.set_aspect('equal', 'box')
        plt.show()

    def find_division_level(self, max_allowed_inflated_outlier):
        divisions = 0
        while True:
            divider = RectangleDivider(self.vertices)
            divider.divide(divisions)
            if divider.sagitta_lateral < max_allowed_inflated_outlier:
                # Update the class members with the final values
                self.rectangles = divider.rectangles
                self.circles = divider.circles
                self.outer_circles = divider.outer_circles
                self.radius = divider.radius
                self.sagitta_lateral = divider.sagitta_lateral
                self.sagitta_longitudinal = divider.sagitta_longitudinal
                self.previous_max_dimension = divider.previous_max_dimension
                return divisions
            divisions += 1
    
    def transform_centers_given_rectangle_pose(self, rectangle_pose):
        rotation_matrix = np.array([
            [np.cos(rectangle_pose[2]), -np.sin(rectangle_pose[2])],
            [np.sin(rectangle_pose[2]), np.cos(rectangle_pose[2])]
        ])
        circle_centers = np.array(self.circles)
        transformed_centers = np.dot(circle_centers, rotation_matrix.T) + rectangle_pose[:2]
        
        return transformed_centers

    def plot_rectangle_and_outlier_centers_with_pose(self, rectangle_pose, plt):
                
        ax = plt.gca()
        
        # Transform the original vertices
        rotation_matrix = np.array([
            [np.cos(rectangle_pose[2]), -np.sin(rectangle_pose[2])],
            [np.sin(rectangle_pose[2]), np.cos(rectangle_pose[2])]
        ])
        transformed_vertices = np.dot(self.vertices, rotation_matrix.T) + rectangle_pose[:2]

        # Draw the transformed rectangle
        rect = plt.Polygon(transformed_vertices, edgecolor='blue', facecolor='none', lw=2)
        ax.add_patch(rect)
        
        # Transform and draw the circle centers
        transformed_centers = self.transform_centers_given_rectangle_pose(rectangle_pose)
        for center_x, center_y in transformed_centers:
            ax.add_patch(plt.Circle((center_x, center_y), self.radius, edgecolor='red', facecolor='none', lw=1, linestyle='--'))

        # ax.set_xlim(np.min(transformed_vertices[:, 0]) - 1, np.max(transformed_vertices[:, 0]) + 1)
        # ax.set_ylim(np.min(transformed_vertices[:, 1]) - 1, np.max(transformed_vertices[:, 1]) + 1)
        # ax.set_aspect('equal', 'box')

# Example usage
if __name__ == "__main__":
    vertices = np.array([
        [0, 0],
        [10, 0],
        [10, 5],
        [0, 5]
    ])
    divider = RectangleDivider(vertices)
    max_allowed_inflated_outlier = 0.5
    division_level = divider.find_division_level(max_allowed_inflated_outlier)
    print(f"Division Level: {division_level}")
    print(f"Radius: {divider.radius}")
    print(f"Sagitta Lateral: {divider.sagitta_lateral}")
    print(f"Sagitta Longitudinal: {divider.sagitta_longitudinal}")
    divider.plot_rectangles_and_circles()
