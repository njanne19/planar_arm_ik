import numpy as np
import matplotlib.pyplot as plt

def draw_simple_arm(ax, theta1, theta2, theta3, l1=1, l2=1, l3=1, x_c=None, y_c=None, x_d=None, y_d=None):
    # Convert angles to radians
    theta1_deg = theta1
    theta2_deg = theta2
    theta3_deg = theta3
    
    theta1 = np.deg2rad(theta1)
    theta2 = np.deg2rad(theta2)
    theta3 = np.deg2rad(theta3)
    
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
    ax.legend(loc='lower right')
    
    # Display angles in the top right-hand corner
    angle_str = f"θ1: {theta1_deg}°\nθ2: {theta2_deg}°\nθ3: {theta3_deg}°\nx_c: {x_c}\ny_c: {y_c}"
    ax.text(0.95, 0.95, angle_str, fontsize=12, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)
    
    # Set equal aspect ratio, grid, and show plot
    ax.set_aspect('equal', 'box')
    ax.grid(True)
    
def draw_armlab_arm(ax, theta1, theta2, theta3, l1=1, l2=1, l3=1, l4 = 1, x_c=None, y_c=None, x_d=None, y_d=None):
    # Convert angles to radians
    theta1_deg = theta1
    theta2_deg = theta2
    theta3_deg = theta3
    
    theta1 = np.deg2rad(theta1)
    theta2 = np.deg2rad(theta2)
    theta3 = np.deg2rad(theta3)
    
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
    ax.legend(loc='lower right')
    
    # Display angles in the top right-hand corner
    angle_str = f"θ1: {theta1_deg}°\nθ2: {theta2_deg}°\nθ3: {theta3_deg}°\nx_c: {x_c}\ny_c: {y_c}"
    ax.text(0.95, 0.95, angle_str, fontsize=12, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)
    
    # Set equal aspect ratio, grid, and show plot
    ax.set_aspect('equal', 'box')
    ax.grid(True)
    
    
    
def main(): 
    # Params for arm 
    l0 = 200 / 1000
    l_gap = 50 / 1000
    l1 = 205.73 / 1000
    l2 = 200 / 1000 
    l3 = (65 + 66 + 43.15/2)/1000

    # Goal points 
    goal_state = [250, -25, 25]
    desired_rotation = 0
    x_d = np.sqrt(np.power(goal_state[0], 2) + np.power(goal_state[1], 2))/1000
    y_d = np.sqrt(np.power(goal_state[0], 2) + np.power(goal_state[2], 2))/1000
    x_d = 0.47
    y_d = 0.2
    x_c = x_d - l3*np.cos(desired_rotation)
    y_c = y_d + l3*np.sin(desired_rotation)
    
    # Detect if point unreachable 
    if np.sqrt(np.power(x_c, 2) + np.power(y_c, 2)) > (l1 + l2):
        print("Point unreachable")
        return 
    elif np.sqrt(np.power(x_c, 2) + np.power(y_c, 2)) < np.abs(l1 - l2):
        print("Point unreachable")
        return

    # Manually change theta1 and theta2 values
    # theta1_simple = 86.589  # angle for l1 in degrees
    # theta2_simple = 150.98  # angle for l2 in degrees
    # theta3_simple = (0 - (theta1_simple + theta2_simple))  # angle for l3 in degrees
    # theta2_input = (np.power(x_c, 2) + np.power(y_c, 2) - np.power(l1, 2) - np.power(l2, 2))/(2*l1*l2)
    # theta2_input_clamped = np.clip(theta2_input, -1, 1)
    
    # theta2_simple = np.arccos(theta2_input_clamped)
    # theta1_simple = np.arctan2(y_c, x_c) - np.arctan2(l2*np.sin((theta2_simple)), (l1 + l2)*np.cos(theta2_simple))  # angle for l1 in degrees
    # theta3_simple = (desired_rotation - (theta1_simple + theta2_simple))  # angle for l3 in degrees
    
    # Chat GPT answers
    theta2_simple = np.arccos((np.power(x_c, 2) + np.power(y_c, 2) - np.power(l1, 2) - np.power(l2, 2))/(2*l1*l2))
    theta2_simple = np.abs(theta2_simple) 
    theta1_simple = np.arctan2(y_c, x_c) - np.arctan2(l2*np.sin((theta2_simple)), (l1 + l2)*np.cos(theta2_simple))
    theta3_simple = (desired_rotation - (theta1_simple + theta2_simple)) 
    # Convert to degrees
    theta1_simple = np.rad2deg(theta1_simple)
    theta2_simple = np.rad2deg(theta2_simple)
    theta3_simple = np.rad2deg(theta3_simple)
    

    # Armlab angles 
    # theta1_armlab = 17.45
    # theta2_armlab = 75.02
    # theta3_armlab = -92.47
    theta1_armlab = 45
    theta2_armlab = 45
    theta3_armlab = 45

    # Generate plot for both 
    arm_figure = plt.figure(figsize=(12, 8))
    arm_axes = arm_figure.subplots(1, 2)

    draw_simple_arm(arm_axes[0], theta1_simple, theta2_simple, theta3_simple, l1, l2, l3, x_c=x_c, y_c=y_c, x_d=x_d, y_d=y_d)
    draw_armlab_arm(arm_axes[1], 0, 0, 0, l0, l_gap, l2, l3, x_c=x_c, y_c=y_c, x_d=x_d, y_d=y_d) 

    # Fix axes scale 
    arm_axes[0].set_xlim([-0.25, 0.75]) 
    arm_axes[0].set_ylim([-0.25, 0.75])
    arm_axes[1].set_xlim([-0.25, 0.75])
    arm_axes[1].set_ylim([-0.25, 0.75])


    plt.show() 


if __name__ == "__main__":
    main()