# imports 
import json, math, logging, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import namedtuple

log = logging.getLogger(__name__)

class PlotException(Exception): 
    """ Custom Exception class with message for this module """
    def __init__(self, value): 
        self.value = value 
        super().__init__(value)
    def __repr__(self): 
        return repr(self.value)


class Controller: 
    """ Class that manages the functionality of the plotting script"""

    def __init__(self, config):
        self.config = config 

    @staticmethod 
    def get_project_root() -> Path:
        return Path(__file__).parent

    @staticmethod
    def get_solution_path() -> Path: 
        return Path(__file__).parent.parent.parent / 'output' / 'solution'

    @staticmethod 
    def get_plot_path() -> Path: 
        return Path(__file__).parent.parent.parent / 'output' / 'plots'

    @staticmethod 
    def get_filename_pattern() -> str: 
        return 'yamal-europe-snapshot'

    @staticmethod
    def load_data() -> namedtuple: 
        dir = Controller.get_solution_path() 
        pattern = Controller.get_filename_pattern() 
        # compute the matrix size and pre-allocate matrices
        num_rows = len(list(Path(dir).glob(f'{pattern}*.json')))
        first_file = f'{dir}/{pattern}-{0}.json'
        with open(first_file, 'r') as f:
            data = json.load(f)
        num_cols =  len(data["initial_pipe_pressure"]["1"]["distance"])
        pressure = np.zeros((num_rows, num_cols))
        flux = np.zeros((num_rows, num_cols))
        time_hrs = np.zeros(num_rows)
        # load data 
        for i in range(0, num_rows):
            file_path = f'{dir}/{pattern}-{i}.json'
            with open(file_path, 'r') as f:
                data = json.load(f)
                pressure[i] = data["initial_pipe_pressure"]["1"]["value"]
                flux[i] = data["initial_pipe_flow"]["1"]["value"]
                time_hrs[i] = data["time"]/3600.0
        # return data
        return namedtuple('Data', ['pressure', 'flux', 'time_hrs'])(pressure, flux, time_hrs)

    @staticmethod 
    def set_rc_params():
        matplotlib.rcParams.update({
            'font.family': 'serif', 
            'font.serif': ['Helvetica'], 
            'font.size': 10
        })

    @staticmethod 
    def plot():
        Controller.set_rc_params()
        pressure, flux, time = Controller.load_data()
        fig, (ax1, ax2, ax3) = plt.subplots(figsize=(5, 6), nrows=3)
        ax1.set_title('Yamal Europe Pipeline Network Solution')
        # plot
        pos1 = ax1.imshow(flux, origin="lower", cmap="RdBu")
        pos2 = ax2.imshow(pressure * 1e-6, origin="lower", cmap="RdBu")
        ax3.plot(time, 1e-6 * pressure[:, -1], "or-")
        # labels
        ax1.set_xlabel("Flux (kg/s)")
        ax1.set_ylabel("Time (hrs)")
        ax1.set_xticks(())
        ax1.set_yticks(())
        # ax1.set_ylim([0, 12])

        ax2.set_xlabel("Pressure (MPa)")
        ax2.set_ylabel("Time (hrs)")
        ax2.set_xticks(())
        ax2.set_yticks(())
        # ax1.set_ylim([0, 12])

        ax3.set_ylabel("Pressure (MPa)")
        ax3.set_xlabel("Time (hrs)")
        ax3.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])

        # color bar positioning
        fig.colorbar(pos1, ax=ax1, location='right', anchor=(0, 0.3), shrink=0.7)
        fig.colorbar(pos2, ax=ax2, location='right', anchor=(0, 0.3), shrink=0.7)
        plt.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)
        plt.savefig(f'{Controller.get_plot_path()}/yamal-europe.png', dpi=600)

def main():
    logging.basicConfig(format='%(asctime)s %(levelname)s--: %(message)s', 
                       level=logging.INFO) 
    try: 
        Controller.plot() 
    except PlotException as pe: 
        log.error(pe)


if __name__ == "__main__":
    main()
