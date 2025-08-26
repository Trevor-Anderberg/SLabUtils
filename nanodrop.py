import numpy as np
import matplotlib.pyplot as plt

class nanodropTSV():
    """
    nanodropTSV(path)
    
    A class for parsing, managing, and visualizing tab-delimited Nanodrop spectrophotometer data files.
    This class reads Nanodrop output files in a specific table format, extracts metadata and absorbance data,
    and provides methods for renaming samples, extracting data for specific samples and wavelengths, and plotting spectra.
    path : str
        Path to the Nanodrop TSV file to be loaded.
    
    Attributes
    ----------
    path : str
        The file path to the loaded Nanodrop data file.
    info : dict
        Dictionary containing header information and lists of date-times, usernames, and sample IDs.
    data : dict
        Dictionary containing:
            - "Wavelengths": numpy array of available wavelengths.
            - Each sample ID as a key, with corresponding absorbance data as a numpy array.
    
    Methods
    -------
    __init__(path)
        Initializes the object, determines file format, and loads data.
    rename_samples(new_names)
    get_data(samples=None, wavelengths=None)
    plot(ax=None, samples=None, wavelengths=None)
        If the file format is invalid or required fields are missing.
        If the length of new sample names does not match the number of samples.

    Notes
    -----
    - Only Nanodrop files with a header containing 'Application:' are supported.
    - If a requested wavelength is not found, the closest available value is used with a warning.
    - If a sample is not found in the data, its row in the output will be filled with NaN or a message will be printed during plotting.
    """

    def __init__(self, path):
        """
        Initializes the object with the given file path, reads the first line to determine the file format,
        and loads data accordingly.
        Parameters:
            path (str): The path to the file to be loaded.
        Raises:
            ValueError: If the file format is invalid (i.e., missing 'Application:' in the header).
        """
        
        self.path = path
        self.info = {}
        self.data = {}
        
        with open(self.path, 'r') as f:
            first_line = f.readline()
            self.is_table = "Application:" in first_line
            if not self.is_table:
                # ideally, we would have a method here for the non table file format
                # until its implementated, we raise an error
                raise ValueError("Invalid file format: missing 'Application:' in header")
            else:
                self._load_table(f)


    def _load_table(self,f):
        
        """
        Loads and parses a tab-delimited table from a file-like object, extracting header information,
        wavelengths, and sample data.
        The method expects the file to start with a header section containing key-value pairs separated
        by ':\t', followed by a line with wavelength information (must include 'Date and Time'), and then
        data rows for each sample.
        
        Parameters
        ----------
        f : file-like object
            An open file or file-like object positioned at the start of the table to be loaded.
        
        Raises
        ------
        ValueError
            If the expected 'Date and Time' field is missing from the wavelength line.
        
        Side Effects
        ------------
        Populates the following instance attributes:
            - self.info: Dictionary with header information and lists of date-times, usernames, and sample IDs.
            - self.data: Dictionary with 'Wavelengths' (numpy array) and sample data (numpy arrays keyed by sample ID).
        """
        


        # return to the beginning of the file
        f.seek(0)
        
        # read info from the header
        while True:
            line = f.readline().strip()
            if ":\t" not in line:
                break
            key, value = line.split(":\t", 1)
            self.info[key.strip()] = value.strip()

        if "Date and Time" not in line:
            raise ValueError("Invalid table format: missing 'Date and Time' in wavelenght line")

        # read wavelengths

        wavelengths = line.split("\t")[3:]
        self.data["Wavelengths"] = np.array(wavelengths, dtype=float)
        
        date_times = []
        usernames = []
        sample_ids = []

        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 4:
                continue
            date_times.append(cols[0].strip())
            sample_ids.append(cols[1].strip())
            usernames.append(cols[2].strip())
            self.data[cols[1].strip()] = np.array(cols[3:], dtype=float)

        self.info["DateTimes"] = date_times
        self.info["Usernames"] = usernames
        self.info["Samples"] = sample_ids

    def rename_samples(self, new_names):
        """
        Renames the samples in the dataset using the provided list of new names.
        
        Parameters
        ----------
        new_names : list of str
            A list containing the new names for the samples. The length of this list
            must match the number of existing samples.
        
        Raises
        ------
        ValueError
            If the length of `new_names` does not match the number of samples.
        
        Notes
        -----
        This method updates both the internal data dictionary and the sample names
        stored in the `info` attribute.
        """

        if len(new_names) != len(self.info["Samples"]):
            raise ValueError("Length of new_names must match number of samples")
        for old_name, new_name in zip(self.info["Samples"], new_names):
            self.data[new_name] = self.data.pop(old_name)
        self.info["Samples"] = new_names

    def get_data(self, samples=None, wavelengths=None):
        """
        Retrieve data for specified samples and wavelengths from the dataset.
        
        Parameters
        ----------
        samples : list, tuple, or None, optional
            List or tuple of sample names to retrieve data for. If None, all samples in self.info["Samples"] are used.
        wavelengths : list, tuple, or None, optional
            Specifies which wavelengths to select:
                - If None, all wavelengths are selected.
                - If a boolean mask (list/array of bools) of the same length as self.data["Wavelengths"], it is used directly.
                - If a tuple, it is interpreted as a (min, max) range and selects wavelengths within that range.
                - If a list of wavelength values, selects those wavelengths. If a value is not found, the closest available wavelength is used with a warning.
        
        Returns
        -------
        data : np.ndarray
            Array of shape (number of samples, number of selected wavelengths) containing the data for the specified samples and wavelengths.
        selected_wavelengths : np.ndarray
            Array of the selected wavelengths corresponding to the columns in `data`.
        
        Raises
        ------
        ValueError
            If a requested wavelength is out of the available range.
        
        Notes
        -----
        - If a requested wavelength is not found in the available wavelengths, the closest value is used and a message is printed.
        - If a sample is not found in the data, its row in the output will be filled with NaN.
        """


        if wavelengths is None:
            mask = np.ones_like(self.data["Wavelengths"], dtype=bool)
        elif len(wavelengths) > 0:
            if issubclass(type(wavelengths[0]), bool) and len(wavelengths) == len(self.data["Wavelengths"]):
                mask = wavelengths
                # if the dtype is boolean, treat wavelengths as a mask
            elif type(wavelengths) is tuple:
                mask = (self.data["Wavelengths"] >= wavelengths[0]) & (self.data["Wavelengths"] <= wavelengths[1])
            else:
                mask = np.isin(self.data["Wavelengths"], wavelengths)
                if not np.all(wavelengths == self.data["Wavelengths"][mask]):
                    # there are wavelengths in the mask that are not in the original list
                    # or things are out of order
                    # manual search
                    mask = []
                    for w in wavelengths:

                        if w < self.data["Wavelengths"][0] or w > self.data["Wavelengths"][-1]:
                            raise ValueError(f"Wavelength {w} is out of range ({self.data['Wavelengths'][0]} - {self.data['Wavelengths'][-1]})")
                        idx = np.argmin(np.abs(self.data["Wavelengths"] - w))
                        if w != self.data["Wavelengths"][idx]:
                            print(f"Wavelength {w} not found, using closest value {self.data['Wavelengths'][idx]}")
                        mask.append(idx)  

        if samples is None:
            samples = self.info["Samples"]
        elif not isinstance(samples, (list, tuple)):
            samples = list(samples)

        if issubclass(type(mask[0]), (bool, np.bool_)):
            num_wavelengths = np.sum(mask)
        else:
            num_wavelengths = len(mask)
        # print(type(mask[0]))
        # print(num_wavelengths)
        # print(mask)

        data = np.zeros((len(samples), num_wavelengths))

        for i, sample in enumerate(samples):
            if sample in self.data:
                data[i] = self.data[sample][mask]
            else:
                data[i] = np.nan

        return data, self.data["Wavelengths"][mask]

    def plot(self, ax=None, samples=None, wavelengths=None):
        """
        Plot absorbance spectra for selected samples and wavelengths.
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The axes on which to plot. If None, a new figure and axes are created.
        samples : list of str, optional
            List of sample names to plot. If None, all samples in `self.info["Samples"]` are plotted.
        wavelengths : array-like, optional
            Wavelength range or mask to plot. If None, all wavelengths are plotted.
            If a list or array of length 2, interpreted as [min, max] wavelength range.
            If a boolean mask of the same length as `self.data["Wavelengths"]`, used directly.
        
        Returns
        -------
        ax : matplotlib.axes.Axes
            The axes with the plotted spectra.
        
        Notes
        -----
        Plots the absorbance spectra for the specified samples and wavelength range.
        Adds labels, title, and legend to the plot.
        Prints a message if a sample is not found in the data.
        """

        if samples is None:
            samples = self.info["Samples"]
        if wavelengths is None:
            mask = np.ones_like(self.data["Wavelengths"], dtype=bool)
        elif len(wavelengths) == 2:
            mask = (self.data["Wavelengths"] >= wavelengths[0]) & (self.data["Wavelengths"] <= wavelengths[1])
        elif len(wavelengths) == len(self.data["Wavelengths"]):
            mask = wavelengths
        for sample in samples:
            if sample in self.data:
                if ax is not None:
                    ax.plot(self.data["Wavelengths"][mask], self.data[sample][mask], label=sample)
                else:
                    plt.plot(self.data["Wavelengths"][mask], self.data[sample][mask], label=sample)
                    ax = plt.gca()
            else:
                print(f"Sample {sample} not found in data.")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Absorbance")
        ax.set_title("Nanodrop Absorbance Spectra")
        ax.grid()
        ax.legend()
        return ax



if __name__ == "__main__":
    
    path = r"S:\tanderberg\Widefield_data\20250825\Nanodrop\UV-Vis 8_25_2025 8_22_28 AM_table.tsv"
    ts = nanodropTSV(path)

    new_names = ['N250_1 S1', 'N250_1 S2', 'D1000 S1', 'D1000 S2', 'N250_1 S3', 'D1000 S3', 'N250_2 S1', 'N250_2 S2', 'N250_2 S3']
               #['Sample 1', 'Sample 2', 'S1 1000 1', 'S2 1000 1', '250 S 3 1', '1000 S 3 1', '250_2 S 1', '250_2 S 2', '250_2 S 3']

    ts.rename_samples(new_names)

    set_1 = ['N250_1 S1', 'N250_1 S2', 'N250_1 S3']
    set_2 = ['D1000 S1', 'D1000 S2', 'D1000 S3']
    set_3 = ['N250_2 S1', 'N250_2 S2', 'N250_2 S3']

    w1 = 500
    w2 = 647

    data_1,_ = ts.get_data(samples=set_1, wavelengths=[w1, w2])
    data_2,_ = ts.get_data(samples=set_2, wavelengths=[w1, w2])
    data_3,_ = ts.get_data(samples=set_3, wavelengths=[w1, w2])

    ratio_1 = data_1[:,1] / data_1[:,0]
    ratio_2 = data_2[:,1] / data_2[:,0]
    ratio_3 = data_3[:,1] / data_3[:,0]

    
    # plt.show()

    fig, axs = plt.subplot_mosaic([['s',"a"],["s",'b'],['s','c']],width_ratios=[1,4], layout='constrained', figsize=(9, 6))
    ax_s = axs["s"]
    ax_a = axs["a"]
    ax_b = axs["b"]
    ax_c = axs["c"]

    ax_s.plot(ratio_2, '--',  alpha=0.5, color='C0')
    ax_s.plot(ratio_1, '--', alpha=0.5, color='C1')
    ax_s.plot(ratio_3, '--', alpha=0.5, color='C2')

    # Make markers fully opaque by re-plotting them on top with alpha=1
    ax_s.plot(ratio_2, 'o', color='C0', alpha=1, label="D1000")
    ax_s.plot(ratio_1, 'o', color='C1', alpha=1, label="N250_1")
    ax_s.plot(ratio_3, 'o', color='C2', alpha=1, label="N250_2")

    ax_s.set_xlabel("Wash Round")
    ax_s.set_ylabel("Absorbance 647/532 Ratio")
    ax_s.set_title("")
    ax_s.legend()
    ax_s.set_yscale("log")
    ax_s.grid()



    ts.plot(ax=ax_a, samples=set_2, wavelengths=[400, 800])
    ax_a.set_title("D1000 " + ax_a.get_title())
    ax_a.axvline(x=w1, color='g', linestyle='--', alpha=0.5)
    ax_a.axvline(x=w2, color='r', linestyle='--', alpha=0.5)

    ts.plot(ax=ax_b, samples=set_1, wavelengths=[400, 800])
    ax_b.set_title("N250_1 " + ax_b.get_title())
    ax_b.axvline(x=w1, color='g', linestyle='--', alpha=0.5)
    ax_b.axvline(x=w2, color='r', linestyle='--', alpha=0.5)

    ts.plot(ax=ax_c, samples=set_3, wavelengths=[400, 800])
    ax_c.set_title("N250_2 " + ax_c.get_title())
    ax_c.axvline(x=w1, color='g', linestyle='--', alpha=0.5)
    ax_c.axvline(x=w2, color='r', linestyle='--', alpha=0.5)
    plt.savefig("nanodrop_plot.png", dpi=300)
    plt.show()




