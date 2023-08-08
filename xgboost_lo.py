import xgboost as xgb
import numpy as np
import joblib  # Import joblib for model serialization

def test_module_load():
    """Display message."""
    return 'This is working!'

def train_xgboost_model(X, y, model_path):
    """
    Train an XGBoost model using the given predictor matrix and output vector,
    and save the model to the specified path.
    
    Parameters:
    X (numpy.ndarray): X by N matrix of predictor values.
    y (numpy.ndarray): 1 by N vector of output values.
    model_path (str): Path to save the trained model.
    
    Returns:
    None
    """
    # Convert input data to DMatrix format
    dmatrix = xgb.DMatrix(X, label=y)
    
    # Set hyperparameters
    params = {
        'objective': 'reg:squarederror',  # You can change this based on your problem
        'eval_metric': 'rmse',  # You can change the evaluation metric
        'max_depth': 3,
        'learning_rate': 0.1,
        'n_estimators': 100,
        # Add more hyperparameters as needed
    }
    
    # Train the XGBoost model
    model = xgb.train(params, dmatrix)
    
    # Save the model to disk
    model.save_model(model_path)

def predict_with_xgboost_model(model_path, X):
    """
    Load a trained XGBoost model from disk and make predictions.
    
    Parameters:
    model_path (str): Path to the saved model.
    X (numpy.ndarray): Input data for prediction.
    
    Returns:
    numpy.ndarray: Predicted output values.
    """
    loaded_model = xgb.Booster()
    loaded_model.load_model(model_path)
    
    dmatrix = xgb.DMatrix(X)
    predictions = loaded_model.predict(dmatrix)
    return predictions