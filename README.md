# StructureSelector-standalone-version
StructureSelector standalone version

```
Usage：perl structure_selector.pl
    -i: path for the input files.
    -t: thredsholds for Puechmaille method[0.5]. split by ; if you want multi values, e.g. 0.5;0.6;0.7
    -o: output path.
    -m: input path contains multi datasets?[0].
    -n: prefix name for the output files.
    -a: input format is admixture or Q-matrix[0]? Default is structure.
    -c: apply chooseK algorithm?
    -l: parse logs for fastStructure or ADMIXTURE.
    -p: group option, popmap files, see examples(can also be e.g. 20,20,20).
    -r: remove ghost method, 'avg' or 'med'.
    -h: help.
  ```
## Example run

for STRUCTURE:
```
perl structure_selector.pl -i input_results_folder -c -o .
```

for ADMIXTURE or other Q-matrix:
```
perl structure_selector.pl -i input_results_folder -p popmap.txt -a -c -o .
```

## Online version
https://lmme.ac.cn/StructureSelector

## Contact

[Yulong Li](mailto:liyulong12@mails.ucas.ac.cn)

Please consider to cite our paper:

> Li YL, Liu JX. 2018 StructureSelector: A web‐based software to select and visualize the optimal number of clusters using multiple methods. Molecular ecology resources, 18(1):176-177. [[link]](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12719)
