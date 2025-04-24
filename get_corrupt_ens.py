import os
import subprocess

def test_file_with_ncdump(file_path):
    result = subprocess.run(
        ['ncdump', '-h', file_path],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print(f"Error in file: {file_path}")
        with open("error_files_ncdump.log", "a") as log_file:
            log_file.write(f"{file_path}\n")
        #return False
    return True

def main():
    base_dir = os.path.join(os.environ['E3SM_ROOT'], 'output', 'UQ', 
                            'UQ_20240104_US-SPR_ICB20TRCNPRDCTCBC')

    # List of valid third-level subfolder names
    valid_subfolders = [
        "T0.00", "T0.00CO2",
        "T2.25", "T2.25CO2",
        "T4.50", "T4.50CO2",
        "T6.75", "T6.75CO2",
        "T9.00", "T9.00CO2",
        "TAMB"
    ]

    # Iterate over second-level folders (g00001 to g04000)
    for second_level in range(1, 4001):
        second_level_path = os.path.join(base_dir, f'g{second_level:05d}')

        for third_level in valid_subfolders:
            e3sm_log = os.path.join(second_level_path, f"e3sm_log_{third_level}.txt")
            with open(e3sm_log) as f:
                log = f.readline()

            if ('error' in log.lower()):
                file_path = second_level_path + ' ' + third_level
                with open("error_files_ncdump.log", "a") as log_file:
                    log_file.write(f"{file_path}\n")

        """# Iterate over third-level folders only if they are in the valid list
        for third_level in valid_subfolders:

            third_level_path = os.path.join(second_level_path, third_level)

            # Check each NetCDF file in the valid third-level folder
            for file in os.listdir(third_level_path):
                if file.endswith('.nc') or file.endswith('.nc4'):
                    file_path = os.path.join(third_level_path, file)
                    result = test_file_with_ncdump(file_path)
                    if not result:
                        # Stop processing further if an error is found
                        return
         """

if __name__ == "__main__":
    main()
