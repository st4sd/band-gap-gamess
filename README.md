# Band-Gap Gamess

__Authors__: James L. McDonagh, Michael Johnston and Hsiang Han Hsu

Virtual experiments for calculating the band-gap of small molecules using [DFT](dft) and [Semi-Empirical](semi-empirical) methods.
Calculations are performed using [GAMESS-US](https://www.msg.chem.iastate.edu/gamess/) 

## Quick links

- [Getting started](#getting-started)
- [Experiment Inputs](#experiment-inputs)
- [Experiment Outputs](#experiment-outputs)
- [Handling Errors](#handling-errors)
- [Computational Details](#computational-details)
- [Development](#development)
- [Help and Support](#help-and-support)
- [Contributing](#contributing)
- [License](#license)

## Getting started

1. Get access to an environment hosting a deployment of the Simulation Toolkit for Scientific Discovery ([ST4SD](https://st4sd.github.io/overview)). 
2. Add the experiments to your ST4SD deployment using [these instructions](#adding-the-experiments-to-st4sd) or check the [ST4SD experiment registry](https://registry.st4sd.res.ibm.com) for other pre-configured variants of the experiments the developers have created.
3. See [here](#experiment-inputs) for the specific inputs these experiments require
4. See [here](http://st4sd.github.io/overview/running-workflows-on-openshift) for general information on running and interacting with virtual-experiments



### Adding the experiments to ST4SD

This will add a basic parameterised package for the Semi-Empirical experiment called `band-gap-pm3-gamess-us`

```python
api.api_experiment_push({"base": {"packages": [{"source": {"git": {"location": {"url": "https://github.com/st4sd/band-gap-gamess.git", "tag": "1.0.0"}}}, "config": {"path": "semi-empirical/homo-lumo-dft-semi-empirical.yaml", "manifestPath": "semi-empirical/manifest.yaml"}}]}, "metadata": {"package": {"name": "band-gap-pm3-gamess-us", "tags": ["latest", "1.0.0"], "maintainer": "https://github.com/michael-johnston", "description": "Uses the PM3 semi-empirical method to perform the geometry optimization and calculate the band-gap and related properties. The calculation is performed with GAMESS-US", "keywords": ["smiles", "computational chemistry", "semi-empirical", "geometry-optimization", "pm3", "homo", "lumo", "band-gap", "gamess-us"]}}, "parameterisation": {"presets": {"runtime": {"args": ["--failSafeDelays=no", "--registerWorkflow=yes"]}}, "executionOptions": {"variables": [{"name": "numberMolecules"}, {"name": "startIndex"}, {"name": "gamess-walltime-minutes"}, {"name": "gamess-grace-period-seconds"}, {"name": "number-processors"}], "platform": ["openshift", "openshift-kubeflux"]}}})
```

This will add a basic parameterised package for the DFT  experiment called `band-gap-dft-gamess-us`

```python
api.api_experiment_push({"metadata": {"package": {"description": "Uses the DFT functional and basis set B3LYP/6-31G(d,p) with Grimme et al's D3 correction to perform geometry optimization and HOMO-LUMO band gap calculation", "keywords": ["smiles", "computational chemistry", "homo-lumo", "dft", "kubeflux"], "maintainer": "https://github.com/michael-johnston", "name": "band-gap-dft-gamess-us", "tags": ["latest", "1.0.0"]}}, "base": {"packages": [{"config": {"manifestPath": "dft/manifest.yaml", "path": "dft/homo-lumo-dft.yaml"}, "name": "main", "source": {"git": {"location": {"tag": "1.0.0", "url": "https://github.com/st4sd/band-gap-gamess.git"}}}}]}, "parameterisation": {"executionOptions": {"data": [], "platform": ["openshift", "openshift-kubeflux"], "variables": [{"name": "numberMolecules"}, {"name": "startIndex"}, {"name": "mem"}, {"name": "gamess-walltime-minutes"}, {"name": "gamess-grace-period-seconds"}, {"name": "number-processors"}]}, "presets": {"environmentVariables": [], "runtime": {"args": ["--failSafeDelays=no", "--registerWorkflow=yes"]}, "variables": [{"name": "functional", "value": "B3LYP"}, {"name": "basis", "value": "GBASIS=N31 NGAUSS=6 NDFUNC=2 NPFUNC=1 DIFFSP=.TRUE. DIFFS=.TRUE."}]}}})
```

You can create your own parameterised packages for these experiments. See [here](https://st4sd.github.io/overview/creating-a-parameterised-package) for more information about how to do this. This repo contains files with the above json in easy to edit form you can use as a basis. 

## Experiment inputs

Both experiments require an input CSV file, called `input_smiles.csv`, with columns `label` and `smiles`. The label column should contain unique integers (e.g. the row number). The `smiles` column should contain the [SMILE](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) representation of the input molecules. The file can contain other columns - these will be ignored. 

Example:

```
label,smiles
0,O=S(=O)([O-])c1c(C(F)(F)F)cc(C(F)(F)F)cc1C(F)(F)F.Cc1cc(OC(C)(C)C)cc(C)c1[S+](c1ccccc1)c1ccccc1
1,O=S(=O)([O-])c1ccc(C(F)(F)F)cc1C(F)(F)F.Cc1ccc([S+](c2ccccc2)c2ccccc2)c(C)c1
```

An example payload for both experiments can be created as follows:

```python
smiles_data = pd.read_csv('smiles_data.csv')[['label', 'SMILES']]
payload = {
    "inputs": [{
       "content": smiles_data.to_csv(index=False),
       "filename": "input_smiles.csv"
   }]
}
```

**Note**: When you submit the `input_smiles.csv` content directly in the payload as in this example, that content can come from a file with any name - here it was `smiles_data.csv`.

By default, the experiment will just measure the first molecule in this file. You can specify a range of molecules to measure using the `startIndex` and `numberMolecules` variables. 

```python
payload = {
    "variables": {
        "startIndex": 0,    #Start at the first molecule in the input file
        "numberMolecules": 2, #Choose two molecules to measure
    },
    "inputs": [{
       "content": smiles_data.to_csv(index=False),
       "filename": "input_smiles.csv"
   }]
}
```

See [here](http://st4sd.github.io/overview/running-workflows-on-openshift#running-a-virtual-experiment) for general information on running virtual-experiments.

### Experiment Outputs

Both experiments have an [interface](http://st4sd.github.io/overview/using-a-virtual-experiment-interface) which allows easy retrieval of the properties they measure. You can see the interface in your ST4SD registry UI when you add one of the experiments (or if you examine other pre-configured packages based on these experiments in the  [global ST4SD experiment registry](https://registry.st4sd.res.ibm.com)) 

See [here](http://st4sd.github.io/overview/running-workflows-on-openshift#retrieving-the-outputs-of-a-virtual-experiment-instance) for more information on accessing the outputs of an experiment run

#### OptimisationResults

The experiment has one key output called `OptimisationResults` 

`OptimisationResults` provides the various HOMO,LUMO and band-gap energies of the molecules. 

To retrieve

```python
path, contents = api.api_rest_uid_output(rest_uid, "OptimisationResults")
print(contents.decode('utf-8'))
```

This would print:

```csv
label,completed,total-energy,homo,lumo,gap,electric-moments,total-time,total-time-per-core
ACOJBQXQMYVURV-UHFFFAOYSA-N_pruned_1,OK,-138248.79125541216,-7.042,-0.724,6.318,5.553946,1129.3,141.16
ACOJBQXQMYVURV-UHFFFAOYSA-N_pruned_1,OK,-138248.52052460602,-6.991,-0.631,6.359,5.496964,1141.9,142.74
```

If the calculation failed for any of the molecules its entries for energetic quantities in this file with be populated with `N/A` (see [Errors](#handling-errors) for more details)


### Handling Errors

There are various errors that can occur during the calculation that prevent the experiment from completing.
The common places to check should an error occur, or you suspect an error has occurred are:

1. `stage0.SMILESToGAMESSInput`
   - Problem: Workflow shutdown on `SMILESToGAMESSInput` component

   - Symptoms:
     - The experiment will finish extremely quickly, typically with an error state referring to a known issue in `SMILESToGAMESSInput`
       For example: `1 jobs failed unexpectedly. Job: stage0.SMILESToGAMESSInput1. Returncode 1. Reason KnownIssue`
     - Can also occur with an unknown error in `SMILESToGAMESSInput`

   - Things to check
     - Check the log file `out.stdout` in `stage1.SMILESToGAMESSInput$N` (see [Computational Details)](#computational-details)). The most likely cause will be a mistake in the SMILES input string. Under the hood RDKit is used to generate 3D structures from each SMILES string followed by a force field optimization.
       if either of these fail this error can occur. The most common error here is a typo in the SMILES string. We suggest that 
        the user tries to visualize the SMILES string for example using [smarts.plus](https://smarts.plus/).

1. `stage1.GeometryOptimisation`

   - Problem: Inability to optimize structure

   - Symptoms: 

     - The experiment will take considerably longer than expected.  
     - The experiment reports success but has a `stage-state` of `component_shutdown`
     - The workflow ends successfully but some homo-lumo energies are missing from the `OptimisationResults` or properties table (instead there appear `N\A`). .
     - The log file reports restart failed - where possible in this case a reason will be given

   - Explanation and things to check
     - Once the GAMESS component (`GeometryOptimisation`) finishes its output is checked by the `AnalyzeEnergies` component. This decides, based on the output file contents, if a successful run is made.  As a default if there is some data to report the output is marked as successful, the `AnalyzeEnergies` component exits with success and the experiment continues. This can lead to a failure at a later stage related to the optimization itself. A common reason for GAMESS failures can be issues writing to disk. A first check should be to make sure the complete output exists in `out.stdout` and in `molecule.dat` in the `GeometryOptimisation` component for the failed molecule ( `molecule.dat` is used to restart with the latest orbitals and coefficients). 

## Computational Details

The computation uses GAMESS version 01. Briefly the calculation consists of:

1. RDKit parsing SMILES strings and generating 3D geometries (`stage0.SMILESToGAMESSInput`)
1. Setting of appropriate GAMESS variables (`stage0.SetBasis`, `stage0.SetFunctional`)
1. GAMESS geometry optimization (`stage1.GeometryOptimisation`)
1. GAMESS output file parsing (`stage1.AnalyzeEnergies`)
1. Output generated to summarize the energy terms (`stage1.ExtractEnergies`)

The identifiers of the workflow components performing each step are given in brackets. These can be used to the various `api.cdb_*` calls to find out more information about these components

For those with knowledge of GAMESS the options set can be found in [`semi-empirical/data-semi-empirical/input_molecule.txt`](semi-empirical/data-semi-empirical/input_molecule.txt). The walltime for the calculation is 1hour. If it is hit the calculation is restarted if suitable data can be found in the `molecule.dat` file under the working directory of the `GeometryOptimisation` component.

A single conformer of the input SMILES string is selected as the lowest energy conformer from 50 generated structures which are minimized by the UFF force field implemented in RDKit

GAMESS does not indicate success/failure of the calculation via its exit-code. When the calculation is completed its consistency is checked by the `AnalyzeEnergies` component. To check the outputs of this for the first molecule in the list of SMILES submitted use: 

```
api.cdb_get_components_last_stdout(instance_uri=m['instance'], component='AnalyzeEnergies0', stage=1)
```

To check the stdout of the GAMESS geometry optimisation for e.g. the first molecule in the list submitted use:

```
api.cdb_get_components_last_stdout(instance_uri=m['instance'], component='GeometryOptimisation0', stage=1)
```

To check the GAMESS output of the `Nth` molecule substitute `N` for `$N` above. 

## Development

1. Fork this repository. You will find 2 virtual experiments in this repository. One that uses [DFT](dft) methods and a second one that uses [Semi-Empirical methods](semi-empirical). 
2. Modify the code
3. Push your code to your forked GitHub repository. Then follow the getting started instructions above.

**Note**: Remember to update your `parameterised package` payload so that it points to your forked GitHub repository.

## Help and Support

Please feel free to reach out to one of the maintainers listed in the [MAINTAINERS.md](MAINTAINERS.md) page.

## Contributing

We always welcome external contributions. Please see our [guidance](CONTRIBUTING.md) for details on how to do so.

## License

This project is licensed under the Apache 2.0 license. Please [see details here](LICENSE.md).
