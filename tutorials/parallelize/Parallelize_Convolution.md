---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: Python [conda env:clonenv]
  language: python
  name: python3
---

# Parallelizing image convolution

+++

## Learning Goals

By the end of this tutorial, you will be able to:

- Employ three parallelization libraries to speed up a serial process.
- Calculate the speedup of the different approaches shown.
- Evaluate which library is suited to your task.

+++

## Introduction

This notebook shows how to speed up an image convolution task using these three libraries:

* Ray: an open-source unified compute framework that makes it easy to scale AI and Python workloads.
* Multiprocessing: part of the standard library; supports spawning processes using an API similar to the threading module; offers both local and remote concurrency, effectively side-stepping the Global Interpreter Lock by using subprocesses instead of threads.
* Dask: developed to natively scale computational packages like numpy, pandas and scikit-learn, and the surrounding ecosystem, to multi-core machines and distributed clusters when datasets exceed memory.

+++

## Imports

* _multiprocessing.Pool_ for multiprocessing using the standard library
* _time_ for timing the processes
* _dask.distributed.Client_ for making a local Dask cluster
* _numpy_ and _scipy.signal_ for numerical work
* _psutil_ for finding the available processors on your machine
* _ray_ for scaling up Python tasks

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install dask[distributed] numpy ray scipy
```

```{code-cell} ipython3
from multiprocessing import Pool
import time

from dask.distributed import Client
import numpy as np
import psutil
import scipy.signal
import ray
```

## Find the cpus available

Find and print the number of cpus
(taken from https://towardsdatascience.com/10x-faster-parallel-python-without-python-multiprocessing-e5017c93cce1)

```{code-cell} ipython3
num_cpus = psutil.cpu_count(logical=True)
print(num_cpus)
```

## Process serially using a conventional loop

+++

Use `scipy.signal` to convolve two 2-dimensional arrays and return a 5x5 downsampled result.

```{code-cell} ipython3
def fconv(image, random_filter):
    return scipy.signal.convolve2d(image, random_filter)[::5, ::5]
```

```{code-cell} ipython3
filters = [np.random.normal(size=(4, 4)) for _ in range(num_cpus)]
```

Process 100 iterations serially, then extrapolate to num_cpus*100

```{code-cell} ipython3
start = time.time()
num_iter = 100
image = np.zeros((3000, 3000))
for i in range(num_iter):
    result = fconv(image, filters[i % num_cpus])
duration_conv = time.time() - start
print("(scaled) conventional duration for {:d} iterations = {:.1f} seconds"
      .format(num_cpus*num_iter, duration_conv*num_cpus))
```

## Process in parallel using Ray

+++

[Documentation for ray](https://docs.ray.io/en/latest/)

+++

The warning raised by `ray.init` only affects shared object usage, which is not an issue for this tutorial. It may harm performance in other scenarios.

```{code-cell} ipython3
ray.init(num_cpus=num_cpus)
```

Use `scipy.signal` to convolve two 2-dimensional arrays and return a 5x5 downsampled result. To use Ray, we decorate the function that is doing the work.

```{code-cell} ipython3
@ray.remote
def fray(image, random_filter):
    return scipy.signal.convolve2d(image, random_filter)[::5, ::5]
```

In the following loop, `ray.put` places the image into shared memory. The call to `ray.get` retrieves the result.

```{code-cell} ipython3
start = time.time()
image = np.zeros((3000, 3000))
for _ in range(100):
    image_id = ray.put(image)
    ray.get([fray.remote(image_id, filters[i]) for i in range(num_cpus)])
duration_ray = time.time() - start
print("Ray duration = {:.1f}, speedup = {:.2f}"
      .format(duration_ray, duration_conv*num_cpus / duration_ray))
```

```{code-cell} ipython3
ray.shutdown()
```

## Process in parallel using multiprocessing

+++

[Documentation for multiprocessing](https://docs.python.org/3/library/multiprocessing.html)

+++

Use `scipy.signal` to convolve two 2-dimensional arrays and return a 5x5 downsampled result. The call to the function has a slightly different form than that for the serial loop.

```{code-cell} ipython3
# Note: Mac and Windows users may need to copy the contents of this cell into a separate '.py' file
# and then import it in order to use the `fmp` function with `multiprocessing`. This has to do with
# differences in what does / does not get copied into the child processes in different operating systems.
import scipy.signal

def fmp(args):
    image, random_filter = args
    return scipy.signal.convolve2d(image, random_filter)[::5, ::5]
```

Use a multiprocessing pool with the number of cpus we found earlier.

```{code-cell} ipython3
pool = Pool(num_cpus)
```

Using `pool.map` is the closest analog in multiprocessing to the Ray API.

```{code-cell} ipython3
start = time.time()
image = np.zeros((3000, 3000))
for _ in range(100):
    pool.map(fmp, zip(num_cpus * [image], filters))
duration_mp = time.time() - start
print("Multiprocessing duration = {:.1f}, speedup = {:.2f}"
      .format(duration_mp, duration_conv*num_cpus / duration_mp))
```

## Process using Dask

+++

[Documentation for Dask](https://www.dask.org/get-started)

+++

Define a Dask distributed client with number of workers set to the number of cpus we found earlier, and with one thread per worker.

```{code-cell} ipython3
client = Client(n_workers=num_cpus, threads_per_worker=1)
```

```{code-cell} ipython3
print(client)
```

Dask recommends scattering the large inputs across the workers, though this makes little difference in execution time.

```{code-cell} ipython3
start = time.time()
image = np.zeros((3000, 3000))
for _ in range(100):
    for j in range(num_cpus):
        big_future = client.scatter((image, filters[j % num_cpus]))
        future = client.submit(fmp, big_future)
duration_dask = time.time() - start
print("Dask duration = {:.1f}, speedup = {:.2f}"
      .format(duration_dask, duration_conv*num_cpus / duration_dask))
```

```{code-cell} ipython3
client.close()
```

## Conclusions

+++

* Ray is the most effective at speeding up the convolution workload by fully utilizing all available processes
* Multiprocessing is second in effectiveness
* Dask delivers the least speedup; perhaps due to having only six processes on the dask.distributed client

+++

## About this notebook

This notebook was developed by David Shupe (shupe@ipac.caltech.edu) in conjunction with Jessica Krick and the IRSA Science Platform team at IPAC.

+++

## Citations

If you use these software packages in your work, please use the following citations:

* Dask: Dask Development Team (2016). Dask: Library for dynamic task scheduling. URL https://dask.org
* Ray: The Ray Development Team. URL https://docs.ray.io

```{code-cell} ipython3

```
