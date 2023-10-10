import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel, QCheckBox, QMainWindow
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import numpy as np 


class App(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.setStyleSheet("""
            QWidget { background-color: #FFFFFF; }
            QLabel { color: #000000; }
            QCheckBox { color: #000000; }
            QSlider::groove:horizontal {
                border: 1px solid #999999;
                height: 8px;
                background: #f0f0f0;
                margin: 2px 0;
            }
            QSlider::handle:horizontal {
                background: #5c5c5c;
                border: 1px solid #5c5c5c;
                width: 18px;
                margin: -2px 0;
                border-radius: 3px;
            }
        """)
        
        self.left = 100
        self.top = 100
        self.title = 'Planar and Armlab IK GUI'
        self.width = 1280
        self.height = 720
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        widget = QWidget(self)
        self.setCentralWidget(widget)

        layout = QVBoxLayout()

        self.canvas = FigureCanvas(plt.figure())
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        self.ax1 = self.canvas.figure.add_subplot(121)
        self.ax2 = self.canvas.figure.add_subplot(122)
        
        sliders_layout = QHBoxLayout()

        # Slider 1 with label above
        slider1_layout = QVBoxLayout()
        self.slider1_label = QLabel("X Position Desired:")
        self.slider1 = QSlider(Qt.Horizontal, self)
        self.slider1.setRange(0, 100)
        self.slider1_default_value = 50
        self.slider1.setValue(self.slider1_default_value)
        self.slider1_readout = QLabel(f"{self.slider1_default_value / 100 * 0.4:5.2f} m")
        self.slider1.valueChanged.connect(lambda: self.slider1_readout.setText(f"{self.slider1.value() / 100 * 0.4:5.2f} m"))
        
        slider1_layout.addWidget(self.slider1_readout)
        slider1_layout.addWidget(self.slider1_label)
        slider1_layout.addWidget(self.slider1)
        sliders_layout.addLayout(slider1_layout)
        

        # Slider 2 with label above
        slider2_layout = QVBoxLayout()
        self.slider2_label = QLabel("Y Position Desired:")
        self.slider2 = QSlider(Qt.Horizontal, self)
        self.slider2.setRange(0, 100)
        self.slider2_default_value = 50
        self.slider2.setValue(self.slider2_default_value)
        self.slider2_readout = QLabel(f"{self.slider2_default_value / 100 * 0.4:5.2f} m")
        self.slider2.valueChanged.connect(lambda: self.slider2_readout.setText(f"{self.slider2.value() / 100 * 0.4:5.2f} m"))
        
        slider2_layout.addWidget(self.slider2_readout)
        slider2_layout.addWidget(self.slider2_label)
        slider2_layout.addWidget(self.slider2)
        sliders_layout.addLayout(slider2_layout)

        # Slider 3 with label above
        slider3_layout = QVBoxLayout()
        self.slider3_label = QLabel("Rotation Desired:")
        self.slider3 = QSlider(Qt.Horizontal, self)
        self.slider3.setRange(-100, 100)
        self.slider3_default_value = 0
        self.slider3.setValue(self.slider3_default_value)
        self.slider3_readout = QLabel(f"{np.rad2deg(self.slider3_default_value / 100 * np.pi):5.2f} deg")
        self.slider3.valueChanged.connect(lambda: self.slider3_readout.setText(f"{np.rad2deg(self.slider3.value() / 100 * np.pi):5.2f} deg"))
        
        slider3_layout.addWidget(self.slider3_readout)
        slider3_layout.addWidget(self.slider3_label)
        slider3_layout.addWidget(self.slider3)
        sliders_layout.addLayout(slider3_layout)

        layout.addLayout(sliders_layout)

        checkbox_layout = QHBoxLayout()
        self.check1 = QCheckBox("Simple Elbow Down", self)
        self.check2 = QCheckBox("Armlab Elbow Down", self)
        checkbox_layout.addWidget(self.check1)
        checkbox_layout.addWidget(self.check2)
        layout.addLayout(checkbox_layout)

        widget.setLayout(layout)

        # Connect events
        self.slider1.valueChanged.connect(self.update_plot)
        self.slider2.valueChanged.connect(self.update_plot)
        self.slider3.valueChanged.connect(self.update_plot)
        self.check1.stateChanged.connect(self.update_plot)
        self.check2.stateChanged.connect(self.update_plot)

        self.update_plot()
        self.show()

    def update_plot(self):
        x_d = (self.slider1.value() / 100) * 0.4
        y_d = (self.slider2.value() / 100) * 0.4
        desired_rotation = (self.slider3.value() / 100) * np.pi
        check1 = self.check1.isChecked()
        check2 = self.check2.isChecked()

        # You can adjust this logic based on your plotting needs
        self.ax1.clear() 
        self.ax2.clear() 
        
        # Base params for arm 
        l0 = 200 / 1000
        l_gap = 50 / 1000
        l1 = 205.73 / 1000
        l2 = 200 / 1000 
        l3 = (65 + 66 + 43.15/2)/1000
        
        # Compute IK params for both arms 
        is_valid_simple, theta1_simple, theta2_simple, theta3_simple, x_c_simple, y_c_simple = recompute_ik_simple(x_d, y_d, desired_rotation, l1, l2, l3, elbow_up=check1)
        is_valid_armlab, theta1_armlab, theta2_armlab, theta3_armlab, x_c_armlab, y_c_armlab = recompute_ik_armlab(x_d, y_d, desired_rotation, l1, l2, l3, elbow_up=check2)
        
        if not is_valid_simple or not is_valid_armlab: 
            mark_unreachable(self.ax1, is_valid_simple, is_valid_armlab) 
            mark_unreachable(self.ax2, is_valid_simple, is_valid_armlab)
        else:           
            # Draw the arm plots 
            draw_simple_arm(self.ax1, theta1_simple, theta2_simple, theta3_simple, l1, l2, l3, x_c=x_c_simple, y_c=y_c_simple, x_d=x_d, y_d=y_d)
            draw_armlab_arm(self.ax2, theta1_armlab, theta2_armlab, theta3_armlab, l0, l_gap, l2, l3, x_c=x_c_armlab, y_c=y_c_armlab, x_d=x_d, y_d=y_d)
            
            # Fix axes scale 
            self.ax1.set_xlim([-0.25, 0.75]) 
            self.ax1.set_ylim([-0.25, 0.75])
            self.ax2.set_xlim([-0.25, 0.75])
            self.ax2.set_ylim([-0.25, 0.75])

        self.canvas.draw()
        
def recompute_ik_simple(x_d, y_d, desired_rotation, l1=1, l2=1, l3=1, elbow_up=True): 

    # Goal points 
    x_c = x_d - l3*np.cos(-1 * desired_rotation)
    y_c = y_d + l3*np.sin(-1 * desired_rotation)
    
    # Detect if point unreachable 
    if np.sqrt(np.power(x_c, 2) + np.power(y_c, 2)) > (l1 + l2):
        return False, None, None, None, None, None
    elif np.sqrt(np.power(x_c, 2) + np.power(y_c, 2)) < np.abs(l1 - l2):
        return False, None, None, None, None, None
    
    #  Two options for theta 2 given elbow up elbow down position
    theta2_simple_1 = np.arccos((np.power(x_c, 2) + np.power(y_c, 2) - np.power(l1, 2) - np.power(l2, 2))/(2*l1*l2))
    theta2_simple_2 = -1 * np.arccos((np.power(x_c, 2) + np.power(y_c, 2) - np.power(l1, 2) - np.power(l2, 2))/(2*l1*l2))
    
    # print(f"theta 2 1: {np.rad2deg(theta2_simple_1)}")
    # print(f"theta 2 2: {np.rad2deg(theta2_simple_2)}")
    
    # Pick which one to use
    if elbow_up: 
        theta2_simple = theta2_simple_1
    else: 
        theta2_simple = theta2_simple_2
    
    theta1_simple = np.arctan(y_c/x_c) - np.arctan((l2*np.sin((theta2_simple)))/(l1 + l2*np.cos(theta2_simple)))
    theta3_simple = (desired_rotation - (theta1_simple + theta2_simple)) 
    
    return True, theta1_simple, theta2_simple, theta3_simple, x_c, y_c

def recompute_ik_armlab(x_d, y_d, desired_rotation, l1=1, l2=1, l3=1, elbow_up=True): 

    # Goal points 
    x_c = x_d - l3*np.cos(-1 * desired_rotation)
    y_c = y_d + l3*np.sin(-1 * desired_rotation)
    
    # Detect if point unreachable 
    if np.sqrt(np.power(x_c, 2) + np.power(y_c, 2)) > (l1 + l2):
        return False, None, None, None, None, None
    elif np.sqrt(np.power(x_c, 2) + np.power(y_c, 2)) < np.abs(l1 - l2):
        return False, None, None, None, None, None
    
    #  Two options for theta 2 given elbow up elbow down position
    theta2_simple_1 = np.arccos((np.power(x_c, 2) + np.power(y_c, 2) - np.power(l1, 2) - np.power(l2, 2))/(2*l1*l2))
    theta2_simple_2 = -1 * np.arccos((np.power(x_c, 2) + np.power(y_c, 2) - np.power(l1, 2) - np.power(l2, 2))/(2*l1*l2))
    
    # print(f"theta 2 1: {np.rad2deg(theta2_simple_1)}")
    # print(f"theta 2 2: {np.rad2deg(theta2_simple_2)}")
    
    # Pick which one to use
    if elbow_up: 
        theta2_simple = theta2_simple_1
    else: 
        theta2_simple = theta2_simple_2
    
    theta1_simple = np.arctan(y_c/x_c) - np.arctan((l2*np.sin((theta2_simple)))/(l1 + l2*np.cos(theta2_simple)))
    theta3_simple = (desired_rotation - (theta1_simple + theta2_simple)) 
    
    # Finally, do correction for armlab arm
    psi = np.arctan(50/200)
    
    # Correction for theta 1
    theta1_corrected = np.pi/2 - psi - theta1_simple 
    
    # Correction for theta 2 
    theta2_corrected = -1*(theta2_simple + (np.pi/2 - psi))
    
    # Correction for theta 3
    theta3_corrected = -1 * theta3_simple
    
    # Finally, check once more if there are any invalid angles
    if np.abs(theta2_corrected) > np.deg2rad(135): 
        return False, None, None, None, None, None
    elif np.abs(theta3_corrected) > np.pi/2: 
        return False, None, None, None, None, None
    
    return True, theta1_corrected, theta2_corrected, theta3_corrected, x_c, y_c

def draw_simple_arm(ax, theta1, theta2, theta3, l1=1, l2=1, l3=1, x_c=None, y_c=None, x_d=None, y_d=None):
    # Convert angles to radians
    # Convert angles to degrees
    theta1_deg = np.rad2deg(theta1)
    theta2_deg = np.rad2deg(theta2)
    theta3_deg = np.rad2deg(theta3)
    # print(f"theta 1: {theta1_deg}")
    # print(f"theta 2: {theta2_deg}")
    # print(f"theta 3: {theta3_deg}")
    
    # Calculate (x, y) positions for first joint
    x1 = l1 * np.cos(theta1)
    y1 = l1 * np.sin(theta1)
    
    # Calculate (x, y) positions for second joint
    x2 = x1 + l2 * np.cos(theta1 + theta2)
    y2 = y1 + l2 * np.sin(theta1 + theta2)
    
    # Calculate (x, y) positions for third joint (end-effector)
    x3 = x2 + l3 * np.cos(theta1 + theta2 + theta3)
    y3 = y2 + l3 * np.sin(theta1 + theta2 + theta3)
    
    # Plot robotic arm
    ax.plot([0, x1], [0, y1], 'o-')
    ax.plot([x1, x2], [y1, y2], 'o-')
    ax.plot([x2, x3], [y2, y3], 'o-')
    
    if x_c is not None and y_c is not None:
        ax.scatter(x_c, y_c, marker='x', color='red', s=100, label='wrist position')
        
    if x_d is not None and y_d is not None:
        ax.scatter(x_d, y_d, marker='o', color='green', s=100, label='desired position')
        
    # Add legend for only the scatter plots
    if x_c is not None or x_d is not None: 
        ax.legend(loc='lower right')
    
    # Display angles in the top right-hand corner
    angle_str = f"θ1: {theta1_deg:5.2f}°\nθ2: {theta2_deg:5.2f}°\nθ3: {theta3_deg:5.2f}°"
    if x_c is not None: 
        angle_str += f"\nx_c: {x_c:5.2f}\ny_c: {y_c:5.2f}"
        
    ax.text(0.95, 0.95, angle_str, fontsize=12, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)
    
    # Set equal aspect ratio, grid, and show plot
    ax.set_aspect('equal', 'box')
    ax.grid(True)
    
    
def draw_armlab_arm(ax, theta1, theta2, theta3, l1=1, l2=1, l3=1, l4 = 1, x_c=None, y_c=None, x_d=None, y_d=None):
    # Convert angles to degrees
    theta1_deg = np.rad2deg(theta1)
    theta2_deg = np.rad2deg(theta2)
    theta3_deg = np.rad2deg(theta3)
    
    # Calculate (x, y) positions for first top node
    x1 = l1 * np.sin(theta1)
    y1 = l1 * np.cos(theta1)
    
    # Calculate (x, y) positions for first joint
    x2 = x1 + l2 * np.cos(theta1)
    y2 = y1 - l2 * np.sin(theta1)
    
    # Calculate position for second joint (position of wrist)
    x3 = x2 + l3 * np.cos(theta1 + theta2)
    y3 = y2 - l3 * np.sin(theta1 + theta2)
    
    # Calculate (x, y) positions for third joint (end-effector)
    x4 = x3 + l4 * np.cos(theta1 + theta2 + theta3)
    y4 = y3 - l4 * np.sin(theta1 + theta2 + theta3)
    
    # Plot robotic arm
    ax.plot([0, x1], [0, y1], 'o-')
    ax.plot([x1, x2], [y1, y2], 'o-') 
    ax.plot([0, x2], [0, y2], 'o-') 
    ax.plot([x2, x3], [y2, y3], 'o-')
    ax.plot([x3, x4], [y3, y4], 'o-')
    
    if x_c is not None and y_c is not None:
        ax.scatter(x_c, y_c, marker='x', color='red', s=100, label='wrist position')
        
    if x_d is not None and y_d is not None:
        ax.scatter(x_d, y_d, marker='o', color='green', s=100, label='desired position')
        
    # Add legend for only the scatter plots
    if x_c is not None or x_d is not None: 
        ax.legend(loc='lower right')
    
    # Display angles in the top right-hand corner
    angle_str = f"θ1: {theta1_deg:5.2f}°\nθ2: {theta2_deg:5.2f}°\nθ3: {theta3_deg:5.2f}°"
    if x_c is not None: 
        angle_str += f"\nx_c: {x_c:5.2f}\ny_c: {y_c:5.2f}"
        
    ax.text(0.95, 0.95, angle_str, fontsize=12, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)
    
    # Set equal aspect ratio, grid, and show plot
    ax.set_aspect('equal', 'box')
    ax.grid(True)
    
def mark_unreachable(ax, reachable_simple, reachable_armlab):
    """
    Draw a red filled square with the text "POINT UNREACHABLE" in the center of the given axis.

    Parameters:
    - ax: A Matplotlib axis object
    """
    
    # Set axis limits
    ax.set_xlim([-0.25, 0.75]) 
    ax.set_ylim([-0.25, 0.75]) 
    
    # Get axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Calculate center of the axis
    center_x = (xlim[0] + xlim[1]) / 2
    center_y = (ylim[0] + ylim[1]) / 2

    # Define the half-length of the square's side
    # (Here, we're taking 10% of the x-axis range, but you can adjust this value)
    half_length = (xlim[1] - xlim[0])

    # Draw the square
    square = plt.Rectangle((center_x - half_length, center_y - half_length), 
                            2*half_length, 2*half_length, color='red')
    ax.add_patch(square)
    
    # Add the text
    if reachable_simple: 
        ax.text(center_x, center_y, 'POINT UNREACHABLE FOR \nARMLAB MODEL ONLY', color='white', 
            weight='bold', ha='center', va='center')
    elif reachable_armlab: 
        ax.text(center_x, center_y, 'POINT UNREACHABLE FOR \nSIMPLE MODEL ONLY', color='white', 
            weight='bold', ha='center', va='center')
    else: 
        ax.text(center_x, center_y, 'POINT UNREACHABLE FOR \nBOTH MODELS', color='white', 
            weight='bold', ha='center', va='center')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
