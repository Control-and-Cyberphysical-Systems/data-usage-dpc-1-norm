# data-usage-dpc-1-norm
Code for generating figures in the paper 'On Data Usage and Predictive Behavior of Data-Driven Predictive Control with 1-Norm Regularization,' published in the IEEE Control Systems Letters Volume 9, pages 943-948. Demonstrates emerging structures related to data usage and predictive behavior in DPC.

*Dependencies:* [MPT 3](https://www.mpt3.org/) is required for parametric optimization, and visualization of piecewise affine functions and polytopes.

## License
The source code in this repository is shared under the MIT license. See the [LICENSE](./LICENSE) file.

## Citing
If you would like to refer to this project in scientific publications, please
reference the original [journal article](https://doi.org/10.1109/LCSYS.2025.3575436)
as:

```bibtex
@ARTICLE{Klaedtke2025dataUsageDPC1Norm,
  author={Kl\"{a}dtke, Manuel and Schulze Darup, Moritz},
  journal={IEEE Control Systems Letters}, 
  title={On Data Usage and Predictive Behavior of Data-Driven Predictive Control With 1-Norm Regularization}, 
  year={2025},
  volume={9},
  number={},
  pages={943-948},
  keywords={Trajectory;Predictive control;Linear systems;Optimization;Predictive models;Noise;Image representation;Costs;Convex hulls;Training;Data-driven control;predictive control;parametric optimization;explainability},
  doi={10.1109/LCSYS.2025.3575436}
}
```


## How to use
Run `GenerateFigures.m`. Doing so will use the same data samples as in the paper. To generate similar figures with new datasets, (un)comment the lines indicated in `GenerateFigures.m`.
