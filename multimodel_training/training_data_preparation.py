from sklearn.model_selection import train_test_split
import os
import polars as pl


class TrainingDataGenerator:
    folder_path = '../data/reference_samples'
    reference_files = [f for f in os.listdir(folder_path) if f.endswith('extended_encoded_reference.csv')]

    for ref_file in reference_files:
        # Data import
        file_path = os.path.join(folder_path, ref_file)
        dtypes = {
            "k-seed": pl.Utf8,
            "rc_seeds": pl.Utf8,
            "runlength_output": pl.Utf8,
            "runlength_rc_output": pl.Utf8,
            "runlength_rc_character": pl.Utf8,
            "runlength_rc_value": pl.Utf8,
            "org_location": pl.Float64,
        }

        # Read the CSV file with explicit column types
        try:
            df = pl.read_csv(file_path, dtypes=dtypes)
        except pl.exceptions.ComputeError as e:
            print(f"Error reading {ref_file}: {e}")
            continue

        # Data cleaning and preparation
        df = df.drop_nulls()  # Drop rows with null values
        df = df.filter(pl.col('org_location') != 0.0)  # Handle rows with 'org_location' == 0.0

        # Splitting features and target
        x = df.select([col for col in df.columns if col != 'org_location'])
        y = df.select(['org_location']).to_series()

        # Splitting data into training, test, and validation sets
        X_train, X_test, y_train, y_test = train_test_split(
            x.to_pandas(), y.to_pandas(), test_size=0.2, shuffle=True, random_state=42
        )
        X_train, X_val, y_train, y_val = train_test_split(
            X_train, y_train, test_size=0.2, shuffle=True, random_state=42
        )

        # Convert split data back to Polars for saving
        X_train = pl.DataFrame(X_train)
        X_test = pl.DataFrame(X_test)
        y_train = pl.DataFrame({'org_location': y_train})
        y_test = pl.DataFrame({'org_location': y_test})
        X_val = pl.DataFrame(X_val)
        y_val = pl.DataFrame({'org_location': y_val})

        # Saving the datasets in the same folder
        base_filename = ref_file.split('_reference.csv')[0]
        X_train.write_csv(os.path.join(folder_path, f"{base_filename}_x_train.csv"))
        X_test.write_csv(os.path.join(folder_path, f"{base_filename}_X_test.csv"))
        y_train.write_csv(os.path.join(folder_path, f"{base_filename}_y_train.csv"))
        y_test.write_csv(os.path.join(folder_path, f"{base_filename}_y_test.csv"))
        X_val.write_csv(os.path.join(folder_path, f"{base_filename}_X_val.csv"))
        y_val.write_csv(os.path.join(folder_path, f"{base_filename}_y_val.csv"))

        # Summary information
        print(f"Processed {ref_file}:")
        print(f"  X training dataset length: {len(X_train)}")
        print(f"  y training dataset length: {len(y_train)}")
