import joblib
import numpy as np


def save_to_bin(data, filename):
    with open(filename, "wb") as file:
        if isinstance(data, np.ndarray):
            np.save(file, data)
        elif isinstance(data, list):
            for array in data:
                np.save(file, array)
        else:
            # Assuming data is an object with NumPy array attributes
            for attr in dir(data):
                attr_value = getattr(data, attr)
                if isinstance(attr_value, np.ndarray):
                    np.save(file, attr_value)


def convert_joblib_to_bin(joblib_filename, bin_filename):
    # Load the model or weights
    data = joblib.load(joblib_filename)

    # Save data to .bin file
    save_to_bin(data, bin_filename)


# Example usage
joblib_filename = "model_weights/centroid_model.joblib"  # Path to your joblib file
bin_filename = "model_weights.bin"  # Desired output .bin file path

convert_joblib_to_bin(joblib_filename, bin_filename)
