

from scipy.linalg import sqrtm
import numpy as np

def simulate_data(n, p, beta_star, sigma=1):
    X = np.random.normal(0, 1, size=(n, p))
    # Generate correlated features
    U = np.random.uniform(size=(p, p))
    X  = X @ U
    eps = np.random.normal(loc=0, scale=sigma, size=n)
    Y = X @ beta_star + eps
    return X, Y


def compute_OLS(X, Y):
    """
    Ordinary Least Squares regression
    """
    # Add a column of ones to X for the intercept
    X = np.column_stack((np.ones(X.shape[0]), X))
    
    # Calculate the coefficients using the Normal Equation
    beta = np.linalg.inv(X.T @ X) @ X.T @ Y
    
    return beta

def compute_combined_OLS(X, Y, X2, Y2, lambda_=0):
    """
    Ordinary Least Squares regression on two datasets
    """
    # Add a column of ones to X for the intercept
    X = np.column_stack((np.ones(X.shape[0]), X))
    X2 = np.column_stack((np.ones(X2.shape[0]), X2))
    # Calculate the coefficients using the Normal Equation
    beta = np.linalg.inv(X.T @ X + lambda_ * X2.T @ X2) @ ( X.T @ Y + lambda_ * X2.T @ Y2)
    
    return beta


X  = np.random.normal(size=(100, 10))






