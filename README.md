# NCBI-Entrez-Query-er
This class object is a simple sequence (RefSeq preferred) 
extractor for the site eutils.ncbi.nlm.nih.gov, which is utility 
development page to create RESTful APIs from these databases.

In order to maximize your experience with the current implementation,
It is strongly suggested that you navigate to 
[https://support.nlm.nih.gov/kbArticle/?pn=KA-05317] and find the link to 
creating an account within the NCBI, so that you can create a private
key. This way the NCBI recognizes your client ip address and will not
throttle you as hard if you make multiple requests to quickly. 

Once you have an API key, you will create a config.json file and input the
following JavaScript notation:

{
	"NCBI_API_KEY":<Your_API_KEY>"
}

From their open NCBIEntrezToFasta.py. You will see more in depth comments
about the function that is currently in the program. File Search "CHANGE THIS"
and follow the instructions of the comment or the layout of the input examples.

In addition, you will need to install biopython, such so you can utilize 
Entrez in the Bio package to your python environment. Assuming your using 
a standard python env (sorry Conda users) you can use PIP (python's package manager) 
to install it. If you do not have PIP installed, run the following from the terminal:
	Windows user: cmd / powershell (do a search in your task bar if you are unfamiliar)
	MacOS/Linux: navigate to your shell (unless you are running no GUI / are headless)
	
	python get-pip.py

After running the installation, from the terminal run the following command:
	pip install biopython

From there, your python environment should be configured and good to start running.
The script has additional input parameters which you need which are gene name and 
sequence count. You can interact with this script from an IDE (Visual Studio, Spyder, Wing, etc.) 
or if you are running a headless terminal run python "directory to where script is saved"\NCBIEntrezToFasta.py
and then enter your inputs. There are some safety parameters to ensure your inputs are valid.

Good Luck!
I have some more functions to add at a later date. 

