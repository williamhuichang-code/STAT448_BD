import numpy as np
import matplotlib.pyplot as plt


def visualize_full_transformation(i_hat_new, j_hat_new, v):
    i_hat_new = np.array(i_hat_new)
    j_hat_new = np.array(j_hat_new)
    v = np.array(v)

    # Build transformation matrix
    A = np.column_stack((i_hat_new, j_hat_new))

    # Transform vector
    v_transformed = A @ v

    # Create grid
    x = np.linspace(-5, 5, 20)
    y = np.linspace(-5, 5, 20)
    X, Y = np.meshgrid(x, y)

    # Flatten grid points
    points = np.vstack([X.ravel(), Y.ravel()])

    # Transform grid
    transformed_points = A @ points
    X_t = transformed_points[0, :].reshape(X.shape)
    Y_t = transformed_points[1, :].reshape(Y.shape)

    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot original grid
    ax.plot(X, Y, color='lightgray', linewidth=0.5)
    ax.plot(X.T, Y.T, color='lightgray', linewidth=0.5)

    # Plot transformed grid
    ax.plot(X_t, Y_t, color='cyan', linewidth=0.8)
    ax.plot(X_t.T, Y_t.T, color='cyan', linewidth=0.8)

    # Original basis
    ax.quiver(0, 0, 1, 0, color='green', scale=1, scale_units='xy', angles='xy')
    ax.quiver(0, 0, 0, 1, color='red', scale=1, scale_units='xy', angles='xy')

    # Transformed basis
    ax.quiver(0, 0, i_hat_new[0], i_hat_new[1],
              color='darkgreen', scale=1, scale_units='xy', angles='xy')
    ax.quiver(0, 0, j_hat_new[0], j_hat_new[1],
              color='darkred', scale=1, scale_units='xy', angles='xy')

    # Original vector
    ax.quiver(0, 0, v[0], v[1],
              color='blue', scale=1, scale_units='xy', angles='xy')

    # Transformed vector
    ax.quiver(0, 0, v_transformed[0], v_transformed[1],
              color='yellow', scale=1, scale_units='xy', angles='xy')

    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_aspect('equal')
    ax.set_title("Linear Transformation of Entire Space")
    ax.grid(False)
    plt.show()


# Example
i_hat_new = (1, -1)
j_hat_new = (2, 0)
v = (-1, 2)

visualize_full_transformation(i_hat_new, j_hat_new, v)