#!/bin/bash
# Conda installation
echo "Installing _libgcc_mutex=0.1=main with conda..."
conda install -y _libgcc_mutex=0.1=main || echo "Conda conflict: _libgcc_mutex=0.1=main" >> conda_conflicts.txt
echo "Installing _openmp_mutex=5.1=1_gnu with conda..."
conda install -y _openmp_mutex=5.1=1_gnu || echo "Conda conflict: _openmp_mutex=5.1=1_gnu" >> conda_conflicts.txt
echo "Installing backcall=0.2.0=pyhd3eb1b0_0 with conda..."
conda install -y backcall=0.2.0=pyhd3eb1b0_0 || echo "Conda conflict: backcall=0.2.0=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing blas=1.0=mkl with conda..."
conda install -y blas=1.0=mkl || echo "Conda conflict: blas=1.0=mkl" >> conda_conflicts.txt
echo "Installing bottleneck=1.3.5=py37h7deecbd_0 with conda..."
conda install -y bottleneck=1.3.5=py37h7deecbd_0 || echo "Conda conflict: bottleneck=1.3.5=py37h7deecbd_0" >> conda_conflicts.txt
echo "Installing brotlipy=0.7.0=py37h27cfd23_1003 with conda..."
conda install -y brotlipy=0.7.0=py37h27cfd23_1003 || echo "Conda conflict: brotlipy=0.7.0=py37h27cfd23_1003" >> conda_conflicts.txt
echo "Installing bzip2=1.0.8=h7b6447c_0 with conda..."
conda install -y bzip2=1.0.8=h7b6447c_0 || echo "Conda conflict: bzip2=1.0.8=h7b6447c_0" >> conda_conflicts.txt
echo "Installing ca-certificates=2022.10.11=h06a4308_0 with conda..."
conda install -y ca-certificates=2022.10.11=h06a4308_0 || echo "Conda conflict: ca-certificates=2022.10.11=h06a4308_0" >> conda_conflicts.txt
echo "Installing cairo=1.16.0=h19f5f5c_2 with conda..."
conda install -y cairo=1.16.0=h19f5f5c_2 || echo "Conda conflict: cairo=1.16.0=h19f5f5c_2" >> conda_conflicts.txt
echo "Installing certifi=2022.9.24=py37h06a4308_0 with conda..."
conda install -y certifi=2022.9.24=py37h06a4308_0 || echo "Conda conflict: certifi=2022.9.24=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing cffi=1.15.1=py37h5eee18b_3 with conda..."
conda install -y cffi=1.15.1=py37h5eee18b_3 || echo "Conda conflict: cffi=1.15.1=py37h5eee18b_3" >> conda_conflicts.txt
echo "Installing charset-normalizer=2.0.4=pyhd3eb1b0_0 with conda..."
conda install -y charset-normalizer=2.0.4=pyhd3eb1b0_0 || echo "Conda conflict: charset-normalizer=2.0.4=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing conda-pack=0.6.0=pyhd3eb1b0_0 with conda..."
conda install -y conda-pack=0.6.0=pyhd3eb1b0_0 || echo "Conda conflict: conda-pack=0.6.0=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing cryptography=38.0.1=py37h9ce1e76_0 with conda..."
conda install -y cryptography=38.0.1=py37h9ce1e76_0 || echo "Conda conflict: cryptography=38.0.1=py37h9ce1e76_0" >> conda_conflicts.txt
echo "Installing cuda=11.6.1=0 with conda..."
conda install -y cuda=11.6.1=0 || echo "Conda conflict: cuda=11.6.1=0" >> conda_conflicts.txt
echo "Installing cuda-cccl=11.6.55=hf6102b2_0 with conda..."
conda install -y cuda-cccl=11.6.55=hf6102b2_0 || echo "Conda conflict: cuda-cccl=11.6.55=hf6102b2_0" >> conda_conflicts.txt
echo "Installing cuda-command-line-tools=11.6.2=0 with conda..."
conda install -y cuda-command-line-tools=11.6.2=0 || echo "Conda conflict: cuda-command-line-tools=11.6.2=0" >> conda_conflicts.txt
echo "Installing cuda-compiler=11.6.2=0 with conda..."
conda install -y cuda-compiler=11.6.2=0 || echo "Conda conflict: cuda-compiler=11.6.2=0" >> conda_conflicts.txt
echo "Installing cuda-cudart=11.6.55=he381448_0 with conda..."
conda install -y cuda-cudart=11.6.55=he381448_0 || echo "Conda conflict: cuda-cudart=11.6.55=he381448_0" >> conda_conflicts.txt
echo "Installing cuda-cudart-dev=11.6.55=h42ad0f4_0 with conda..."
conda install -y cuda-cudart-dev=11.6.55=h42ad0f4_0 || echo "Conda conflict: cuda-cudart-dev=11.6.55=h42ad0f4_0" >> conda_conflicts.txt
echo "Installing cuda-cuobjdump=11.6.124=h2eeebcb_0 with conda..."
conda install -y cuda-cuobjdump=11.6.124=h2eeebcb_0 || echo "Conda conflict: cuda-cuobjdump=11.6.124=h2eeebcb_0" >> conda_conflicts.txt
echo "Installing cuda-cupti=11.6.124=h86345e5_0 with conda..."
conda install -y cuda-cupti=11.6.124=h86345e5_0 || echo "Conda conflict: cuda-cupti=11.6.124=h86345e5_0" >> conda_conflicts.txt
echo "Installing cuda-cuxxfilt=11.6.124=hecbf4f6_0 with conda..."
conda install -y cuda-cuxxfilt=11.6.124=hecbf4f6_0 || echo "Conda conflict: cuda-cuxxfilt=11.6.124=hecbf4f6_0" >> conda_conflicts.txt
echo "Installing cuda-driver-dev=11.6.55=0 with conda..."
conda install -y cuda-driver-dev=11.6.55=0 || echo "Conda conflict: cuda-driver-dev=11.6.55=0" >> conda_conflicts.txt
echo "Installing cuda-gdb=12.0.90=0 with conda..."
conda install -y cuda-gdb=12.0.90=0 || echo "Conda conflict: cuda-gdb=12.0.90=0" >> conda_conflicts.txt
echo "Installing cuda-libraries=11.6.1=0 with conda..."
conda install -y cuda-libraries=11.6.1=0 || echo "Conda conflict: cuda-libraries=11.6.1=0" >> conda_conflicts.txt
echo "Installing cuda-libraries-dev=11.6.1=0 with conda..."
conda install -y cuda-libraries-dev=11.6.1=0 || echo "Conda conflict: cuda-libraries-dev=11.6.1=0" >> conda_conflicts.txt
echo "Installing cuda-memcheck=11.8.86=0 with conda..."
conda install -y cuda-memcheck=11.8.86=0 || echo "Conda conflict: cuda-memcheck=11.8.86=0" >> conda_conflicts.txt
echo "Installing cuda-nsight=12.0.78=0 with conda..."
conda install -y cuda-nsight=12.0.78=0 || echo "Conda conflict: cuda-nsight=12.0.78=0" >> conda_conflicts.txt
echo "Installing cuda-nsight-compute=12.0.0=0 with conda..."
conda install -y cuda-nsight-compute=12.0.0=0 || echo "Conda conflict: cuda-nsight-compute=12.0.0=0" >> conda_conflicts.txt
echo "Installing cuda-nvcc=11.6.124=hbba6d2d_0 with conda..."
conda install -y cuda-nvcc=11.6.124=hbba6d2d_0 || echo "Conda conflict: cuda-nvcc=11.6.124=hbba6d2d_0" >> conda_conflicts.txt
echo "Installing cuda-nvdisasm=12.0.76=0 with conda..."
conda install -y cuda-nvdisasm=12.0.76=0 || echo "Conda conflict: cuda-nvdisasm=12.0.76=0" >> conda_conflicts.txt
echo "Installing cuda-nvml-dev=11.6.55=haa9ef22_0 with conda..."
conda install -y cuda-nvml-dev=11.6.55=haa9ef22_0 || echo "Conda conflict: cuda-nvml-dev=11.6.55=haa9ef22_0" >> conda_conflicts.txt
echo "Installing cuda-nvprof=12.0.90=0 with conda..."
conda install -y cuda-nvprof=12.0.90=0 || echo "Conda conflict: cuda-nvprof=12.0.90=0" >> conda_conflicts.txt
echo "Installing cuda-nvprune=11.6.124=he22ec0a_0 with conda..."
conda install -y cuda-nvprune=11.6.124=he22ec0a_0 || echo "Conda conflict: cuda-nvprune=11.6.124=he22ec0a_0" >> conda_conflicts.txt
echo "Installing cuda-nvrtc=11.6.124=h020bade_0 with conda..."
conda install -y cuda-nvrtc=11.6.124=h020bade_0 || echo "Conda conflict: cuda-nvrtc=11.6.124=h020bade_0" >> conda_conflicts.txt
echo "Installing cuda-nvrtc-dev=11.6.124=h249d397_0 with conda..."
conda install -y cuda-nvrtc-dev=11.6.124=h249d397_0 || echo "Conda conflict: cuda-nvrtc-dev=11.6.124=h249d397_0" >> conda_conflicts.txt
echo "Installing cuda-nvtx=11.6.124=h0630a44_0 with conda..."
conda install -y cuda-nvtx=11.6.124=h0630a44_0 || echo "Conda conflict: cuda-nvtx=11.6.124=h0630a44_0" >> conda_conflicts.txt
echo "Installing cuda-nvvp=12.0.90=0 with conda..."
conda install -y cuda-nvvp=12.0.90=0 || echo "Conda conflict: cuda-nvvp=12.0.90=0" >> conda_conflicts.txt
echo "Installing cuda-runtime=11.6.1=0 with conda..."
conda install -y cuda-runtime=11.6.1=0 || echo "Conda conflict: cuda-runtime=11.6.1=0" >> conda_conflicts.txt
echo "Installing cuda-samples=11.6.101=h8efea70_0 with conda..."
conda install -y cuda-samples=11.6.101=h8efea70_0 || echo "Conda conflict: cuda-samples=11.6.101=h8efea70_0" >> conda_conflicts.txt
echo "Installing cuda-sanitizer-api=12.0.90=0 with conda..."
conda install -y cuda-sanitizer-api=12.0.90=0 || echo "Conda conflict: cuda-sanitizer-api=12.0.90=0" >> conda_conflicts.txt
echo "Installing cuda-toolkit=11.6.1=0 with conda..."
conda install -y cuda-toolkit=11.6.1=0 || echo "Conda conflict: cuda-toolkit=11.6.1=0" >> conda_conflicts.txt
echo "Installing cuda-tools=11.6.1=0 with conda..."
conda install -y cuda-tools=11.6.1=0 || echo "Conda conflict: cuda-tools=11.6.1=0" >> conda_conflicts.txt
echo "Installing cuda-visual-tools=11.6.1=0 with conda..."
conda install -y cuda-visual-tools=11.6.1=0 || echo "Conda conflict: cuda-visual-tools=11.6.1=0" >> conda_conflicts.txt
echo "Installing debugpy=1.5.1=py37h295c915_0 with conda..."
conda install -y debugpy=1.5.1=py37h295c915_0 || echo "Conda conflict: debugpy=1.5.1=py37h295c915_0" >> conda_conflicts.txt
echo "Installing decorator=5.1.1=pyhd3eb1b0_0 with conda..."
conda install -y decorator=5.1.1=pyhd3eb1b0_0 || echo "Conda conflict: decorator=5.1.1=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing entrypoints=0.4=py37h06a4308_0 with conda..."
conda install -y entrypoints=0.4=py37h06a4308_0 || echo "Conda conflict: entrypoints=0.4=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing ffmpeg=4.3=hf484d3e_0 with conda..."
conda install -y ffmpeg=4.3=hf484d3e_0 || echo "Conda conflict: ffmpeg=4.3=hf484d3e_0" >> conda_conflicts.txt
echo "Installing flit-core=3.6.0=pyhd3eb1b0_0 with conda..."
conda install -y flit-core=3.6.0=pyhd3eb1b0_0 || echo "Conda conflict: flit-core=3.6.0=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing fontconfig=2.14.1=hef1e5e3_0 with conda..."
conda install -y fontconfig=2.14.1=hef1e5e3_0 || echo "Conda conflict: fontconfig=2.14.1=hef1e5e3_0" >> conda_conflicts.txt
echo "Installing freetype=2.12.1=h4a9f257_0 with conda..."
conda install -y freetype=2.12.1=h4a9f257_0 || echo "Conda conflict: freetype=2.12.1=h4a9f257_0" >> conda_conflicts.txt
echo "Installing gds-tools=1.5.0.59=0 with conda..."
conda install -y gds-tools=1.5.0.59=0 || echo "Conda conflict: gds-tools=1.5.0.59=0" >> conda_conflicts.txt
echo "Installing giflib=5.2.1=h7b6447c_0 with conda..."
conda install -y giflib=5.2.1=h7b6447c_0 || echo "Conda conflict: giflib=5.2.1=h7b6447c_0" >> conda_conflicts.txt
echo "Installing glib=2.69.1=he621ea3_2 with conda..."
conda install -y glib=2.69.1=he621ea3_2 || echo "Conda conflict: glib=2.69.1=he621ea3_2" >> conda_conflicts.txt
echo "Installing gmp=6.2.1=h295c915_3 with conda..."
conda install -y gmp=6.2.1=h295c915_3 || echo "Conda conflict: gmp=6.2.1=h295c915_3" >> conda_conflicts.txt
echo "Installing gnutls=3.6.15=he1e5248_0 with conda..."
conda install -y gnutls=3.6.15=he1e5248_0 || echo "Conda conflict: gnutls=3.6.15=he1e5248_0" >> conda_conflicts.txt
echo "Installing icu=58.2=he6710b0_3 with conda..."
conda install -y icu=58.2=he6710b0_3 || echo "Conda conflict: icu=58.2=he6710b0_3" >> conda_conflicts.txt
echo "Installing idna=3.4=py37h06a4308_0 with conda..."
conda install -y idna=3.4=py37h06a4308_0 || echo "Conda conflict: idna=3.4=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing intel-openmp=2021.4.0=h06a4308_3561 with conda..."
conda install -y intel-openmp=2021.4.0=h06a4308_3561 || echo "Conda conflict: intel-openmp=2021.4.0=h06a4308_3561" >> conda_conflicts.txt
echo "Installing ipykernel=6.15.2=py37h06a4308_0 with conda..."
conda install -y ipykernel=6.15.2=py37h06a4308_0 || echo "Conda conflict: ipykernel=6.15.2=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing jpeg=9e=h7f8727e_0 with conda..."
conda install -y jpeg=9e=h7f8727e_0 || echo "Conda conflict: jpeg=9e=h7f8727e_0" >> conda_conflicts.txt
echo "Installing jupyter_client=7.4.7=py37h06a4308_0 with conda..."
conda install -y jupyter_client=7.4.7=py37h06a4308_0 || echo "Conda conflict: jupyter_client=7.4.7=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing jupyter_core=4.11.2=py37h06a4308_0 with conda..."
conda install -y jupyter_core=4.11.2=py37h06a4308_0 || echo "Conda conflict: jupyter_core=4.11.2=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing lame=3.100=h7b6447c_0 with conda..."
conda install -y lame=3.100=h7b6447c_0 || echo "Conda conflict: lame=3.100=h7b6447c_0" >> conda_conflicts.txt
echo "Installing lcms2=2.12=h3be6417_0 with conda..."
conda install -y lcms2=2.12=h3be6417_0 || echo "Conda conflict: lcms2=2.12=h3be6417_0" >> conda_conflicts.txt
echo "Installing ld_impl_linux-64=2.38=h1181459_1 with conda..."
conda install -y ld_impl_linux-64=2.38=h1181459_1 || echo "Conda conflict: ld_impl_linux-64=2.38=h1181459_1" >> conda_conflicts.txt
echo "Installing lerc=3.0=h295c915_0 with conda..."
conda install -y lerc=3.0=h295c915_0 || echo "Conda conflict: lerc=3.0=h295c915_0" >> conda_conflicts.txt
echo "Installing libboost=1.73.0=h28710b8_12 with conda..."
conda install -y libboost=1.73.0=h28710b8_12 || echo "Conda conflict: libboost=1.73.0=h28710b8_12" >> conda_conflicts.txt
echo "Installing libcublas=11.9.2.110=h5e84587_0 with conda..."
conda install -y libcublas=11.9.2.110=h5e84587_0 || echo "Conda conflict: libcublas=11.9.2.110=h5e84587_0" >> conda_conflicts.txt
echo "Installing libcublas-dev=11.9.2.110=h5c901ab_0 with conda..."
conda install -y libcublas-dev=11.9.2.110=h5c901ab_0 || echo "Conda conflict: libcublas-dev=11.9.2.110=h5c901ab_0" >> conda_conflicts.txt
echo "Installing libcufft=10.7.1.112=hf425ae0_0 with conda..."
conda install -y libcufft=10.7.1.112=hf425ae0_0 || echo "Conda conflict: libcufft=10.7.1.112=hf425ae0_0" >> conda_conflicts.txt
echo "Installing libcufft-dev=10.7.1.112=ha5ce4c0_0 with conda..."
conda install -y libcufft-dev=10.7.1.112=ha5ce4c0_0 || echo "Conda conflict: libcufft-dev=10.7.1.112=ha5ce4c0_0" >> conda_conflicts.txt
echo "Installing libcufile=1.5.0.59=0 with conda..."
conda install -y libcufile=1.5.0.59=0 || echo "Conda conflict: libcufile=1.5.0.59=0" >> conda_conflicts.txt
echo "Installing libcufile-dev=1.5.0.59=0 with conda..."
conda install -y libcufile-dev=1.5.0.59=0 || echo "Conda conflict: libcufile-dev=1.5.0.59=0" >> conda_conflicts.txt
echo "Installing libcurand=10.3.1.50=0 with conda..."
conda install -y libcurand=10.3.1.50=0 || echo "Conda conflict: libcurand=10.3.1.50=0" >> conda_conflicts.txt
echo "Installing libcurand-dev=10.3.1.50=0 with conda..."
conda install -y libcurand-dev=10.3.1.50=0 || echo "Conda conflict: libcurand-dev=10.3.1.50=0" >> conda_conflicts.txt
echo "Installing libcusolver=11.3.4.124=h33c3c4e_0 with conda..."
conda install -y libcusolver=11.3.4.124=h33c3c4e_0 || echo "Conda conflict: libcusolver=11.3.4.124=h33c3c4e_0" >> conda_conflicts.txt
echo "Installing libcusparse=11.7.2.124=h7538f96_0 with conda..."
conda install -y libcusparse=11.7.2.124=h7538f96_0 || echo "Conda conflict: libcusparse=11.7.2.124=h7538f96_0" >> conda_conflicts.txt
echo "Installing libcusparse-dev=11.7.2.124=hbbe9722_0 with conda..."
conda install -y libcusparse-dev=11.7.2.124=hbbe9722_0 || echo "Conda conflict: libcusparse-dev=11.7.2.124=hbbe9722_0" >> conda_conflicts.txt
echo "Installing libdeflate=1.8=h7f8727e_5 with conda..."
conda install -y libdeflate=1.8=h7f8727e_5 || echo "Conda conflict: libdeflate=1.8=h7f8727e_5" >> conda_conflicts.txt
echo "Installing libffi=3.4.2=h6a678d5_6 with conda..."
conda install -y libffi=3.4.2=h6a678d5_6 || echo "Conda conflict: libffi=3.4.2=h6a678d5_6" >> conda_conflicts.txt
echo "Installing libgcc-ng=11.2.0=h1234567_1 with conda..."
conda install -y libgcc-ng=11.2.0=h1234567_1 || echo "Conda conflict: libgcc-ng=11.2.0=h1234567_1" >> conda_conflicts.txt
echo "Installing libgomp=11.2.0=h1234567_1 with conda..."
conda install -y libgomp=11.2.0=h1234567_1 || echo "Conda conflict: libgomp=11.2.0=h1234567_1" >> conda_conflicts.txt
echo "Installing libiconv=1.16=h7f8727e_2 with conda..."
conda install -y libiconv=1.16=h7f8727e_2 || echo "Conda conflict: libiconv=1.16=h7f8727e_2" >> conda_conflicts.txt
echo "Installing libidn2=2.3.2=h7f8727e_0 with conda..."
conda install -y libidn2=2.3.2=h7f8727e_0 || echo "Conda conflict: libidn2=2.3.2=h7f8727e_0" >> conda_conflicts.txt
echo "Installing libnpp=11.6.3.124=hd2722f0_0 with conda..."
conda install -y libnpp=11.6.3.124=hd2722f0_0 || echo "Conda conflict: libnpp=11.6.3.124=hd2722f0_0" >> conda_conflicts.txt
echo "Installing libnpp-dev=11.6.3.124=h3c42840_0 with conda..."
conda install -y libnpp-dev=11.6.3.124=h3c42840_0 || echo "Conda conflict: libnpp-dev=11.6.3.124=h3c42840_0" >> conda_conflicts.txt
echo "Installing libnvjpeg=11.6.2.124=hd473ad6_0 with conda..."
conda install -y libnvjpeg=11.6.2.124=hd473ad6_0 || echo "Conda conflict: libnvjpeg=11.6.2.124=hd473ad6_0" >> conda_conflicts.txt
echo "Installing libnvjpeg-dev=11.6.2.124=hb5906b9_0 with conda..."
conda install -y libnvjpeg-dev=11.6.2.124=hb5906b9_0 || echo "Conda conflict: libnvjpeg-dev=11.6.2.124=hb5906b9_0" >> conda_conflicts.txt
echo "Installing libpng=1.6.37=hbc83047_0 with conda..."
conda install -y libpng=1.6.37=hbc83047_0 || echo "Conda conflict: libpng=1.6.37=hbc83047_0" >> conda_conflicts.txt
echo "Installing libsodium=1.0.18=h7b6447c_0 with conda..."
conda install -y libsodium=1.0.18=h7b6447c_0 || echo "Conda conflict: libsodium=1.0.18=h7b6447c_0" >> conda_conflicts.txt
echo "Installing libstdcxx-ng=11.2.0=h1234567_1 with conda..."
conda install -y libstdcxx-ng=11.2.0=h1234567_1 || echo "Conda conflict: libstdcxx-ng=11.2.0=h1234567_1" >> conda_conflicts.txt
echo "Installing libtasn1=4.16.0=h27cfd23_0 with conda..."
conda install -y libtasn1=4.16.0=h27cfd23_0 || echo "Conda conflict: libtasn1=4.16.0=h27cfd23_0" >> conda_conflicts.txt
echo "Installing libtiff=4.4.0=hecacb30_2 with conda..."
conda install -y libtiff=4.4.0=hecacb30_2 || echo "Conda conflict: libtiff=4.4.0=hecacb30_2" >> conda_conflicts.txt
echo "Installing libunistring=0.9.10=h27cfd23_0 with conda..."
conda install -y libunistring=0.9.10=h27cfd23_0 || echo "Conda conflict: libunistring=0.9.10=h27cfd23_0" >> conda_conflicts.txt
echo "Installing libwebp=1.2.4=h11a3e52_0 with conda..."
conda install -y libwebp=1.2.4=h11a3e52_0 || echo "Conda conflict: libwebp=1.2.4=h11a3e52_0" >> conda_conflicts.txt
echo "Installing libwebp-base=1.2.4=h5eee18b_0 with conda..."
conda install -y libwebp-base=1.2.4=h5eee18b_0 || echo "Conda conflict: libwebp-base=1.2.4=h5eee18b_0" >> conda_conflicts.txt
echo "Installing libxcb=1.15=h7f8727e_0 with conda..."
conda install -y libxcb=1.15=h7f8727e_0 || echo "Conda conflict: libxcb=1.15=h7f8727e_0" >> conda_conflicts.txt
echo "Installing libxml2=2.9.14=h74e7548_0 with conda..."
conda install -y libxml2=2.9.14=h74e7548_0 || echo "Conda conflict: libxml2=2.9.14=h74e7548_0" >> conda_conflicts.txt
echo "Installing lz4-c=1.9.4=h6a678d5_0 with conda..."
conda install -y lz4-c=1.9.4=h6a678d5_0 || echo "Conda conflict: lz4-c=1.9.4=h6a678d5_0" >> conda_conflicts.txt
echo "Installing matplotlib-inline=0.1.6=py37h06a4308_0 with conda..."
conda install -y matplotlib-inline=0.1.6=py37h06a4308_0 || echo "Conda conflict: matplotlib-inline=0.1.6=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing mkl=2021.4.0=h06a4308_640 with conda..."
conda install -y mkl=2021.4.0=h06a4308_640 || echo "Conda conflict: mkl=2021.4.0=h06a4308_640" >> conda_conflicts.txt
echo "Installing mkl-service=2.4.0=py37h7f8727e_0 with conda..."
conda install -y mkl-service=2.4.0=py37h7f8727e_0 || echo "Conda conflict: mkl-service=2.4.0=py37h7f8727e_0" >> conda_conflicts.txt
echo "Installing mkl_fft=1.3.1=py37hd3c417c_0 with conda..."
conda install -y mkl_fft=1.3.1=py37hd3c417c_0 || echo "Conda conflict: mkl_fft=1.3.1=py37hd3c417c_0" >> conda_conflicts.txt
echo "Installing mkl_random=1.2.2=py37h51133e4_0 with conda..."
conda install -y mkl_random=1.2.2=py37h51133e4_0 || echo "Conda conflict: mkl_random=1.2.2=py37h51133e4_0" >> conda_conflicts.txt
echo "Installing ncurses=6.3=h5eee18b_3 with conda..."
conda install -y ncurses=6.3=h5eee18b_3 || echo "Conda conflict: ncurses=6.3=h5eee18b_3" >> conda_conflicts.txt
echo "Installing nest-asyncio=1.5.5=py37h06a4308_0 with conda..."
conda install -y nest-asyncio=1.5.5=py37h06a4308_0 || echo "Conda conflict: nest-asyncio=1.5.5=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing nettle=3.7.3=hbbd107a_1 with conda..."
conda install -y nettle=3.7.3=hbbd107a_1 || echo "Conda conflict: nettle=3.7.3=hbbd107a_1" >> conda_conflicts.txt
echo "Installing nsight-compute=2022.4.0.15=0 with conda..."
conda install -y nsight-compute=2022.4.0.15=0 || echo "Conda conflict: nsight-compute=2022.4.0.15=0" >> conda_conflicts.txt
echo "Installing numexpr=2.8.4=py37he184ba9_0 with conda..."
conda install -y numexpr=2.8.4=py37he184ba9_0 || echo "Conda conflict: numexpr=2.8.4=py37he184ba9_0" >> conda_conflicts.txt
echo "Installing numpy=1.21.5=py37h6c91a56_3 with conda..."
conda install -y numpy=1.21.5=py37h6c91a56_3 || echo "Conda conflict: numpy=1.21.5=py37h6c91a56_3" >> conda_conflicts.txt
echo "Installing numpy-base=1.21.5=py37ha15fc14_3 with conda..."
conda install -y numpy-base=1.21.5=py37ha15fc14_3 || echo "Conda conflict: numpy-base=1.21.5=py37ha15fc14_3" >> conda_conflicts.txt
echo "Installing openh264=2.1.1=h4ff587b_0 with conda..."
conda install -y openh264=2.1.1=h4ff587b_0 || echo "Conda conflict: openh264=2.1.1=h4ff587b_0" >> conda_conflicts.txt
echo "Installing openssl=1.1.1s=h7f8727e_0 with conda..."
conda install -y openssl=1.1.1s=h7f8727e_0 || echo "Conda conflict: openssl=1.1.1s=h7f8727e_0" >> conda_conflicts.txt
echo "Installing packaging=21.3=pyhd3eb1b0_0 with conda..."
conda install -y packaging=21.3=pyhd3eb1b0_0 || echo "Conda conflict: packaging=21.3=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing pandas=1.3.5=py37h8c16a72_0 with conda..."
conda install -y pandas=1.3.5=py37h8c16a72_0 || echo "Conda conflict: pandas=1.3.5=py37h8c16a72_0" >> conda_conflicts.txt
echo "Installing parso=0.8.3=pyhd3eb1b0_0 with conda..."
conda install -y parso=0.8.3=pyhd3eb1b0_0 || echo "Conda conflict: parso=0.8.3=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing pcre=8.45=h295c915_0 with conda..."
conda install -y pcre=8.45=h295c915_0 || echo "Conda conflict: pcre=8.45=h295c915_0" >> conda_conflicts.txt
echo "Installing pexpect=4.8.0=pyhd3eb1b0_3 with conda..."
conda install -y pexpect=4.8.0=pyhd3eb1b0_3 || echo "Conda conflict: pexpect=4.8.0=pyhd3eb1b0_3" >> conda_conflicts.txt
echo "Installing pickleshare=0.7.5=pyhd3eb1b0_1003 with conda..."
conda install -y pickleshare=0.7.5=pyhd3eb1b0_1003 || echo "Conda conflict: pickleshare=0.7.5=pyhd3eb1b0_1003" >> conda_conflicts.txt
echo "Installing pillow=9.2.0=py37hace64e9_1 with conda..."
conda install -y pillow=9.2.0=py37hace64e9_1 || echo "Conda conflict: pillow=9.2.0=py37hace64e9_1" >> conda_conflicts.txt
echo "Installing pip=22.3.1=py37h06a4308_0 with conda..."
conda install -y pip=22.3.1=py37h06a4308_0 || echo "Conda conflict: pip=22.3.1=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing pixman=0.40.0=h7f8727e_1 with conda..."
conda install -y pixman=0.40.0=h7f8727e_1 || echo "Conda conflict: pixman=0.40.0=h7f8727e_1" >> conda_conflicts.txt
echo "Installing ptyprocess=0.7.0=pyhd3eb1b0_2 with conda..."
conda install -y ptyprocess=0.7.0=pyhd3eb1b0_2 || echo "Conda conflict: ptyprocess=0.7.0=pyhd3eb1b0_2" >> conda_conflicts.txt
echo "Installing py-boost=1.73.0=py37h51133e4_12 with conda..."
conda install -y py-boost=1.73.0=py37h51133e4_12 || echo "Conda conflict: py-boost=1.73.0=py37h51133e4_12" >> conda_conflicts.txt
echo "Installing pycparser=2.21=pyhd3eb1b0_0 with conda..."
conda install -y pycparser=2.21=pyhd3eb1b0_0 || echo "Conda conflict: pycparser=2.21=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing pyopenssl=22.0.0=pyhd3eb1b0_0 with conda..."
conda install -y pyopenssl=22.0.0=pyhd3eb1b0_0 || echo "Conda conflict: pyopenssl=22.0.0=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing pyparsing=3.0.9=py37h06a4308_0 with conda..."
conda install -y pyparsing=3.0.9=py37h06a4308_0 || echo "Conda conflict: pyparsing=3.0.9=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing pysocks=1.7.1=py37_1 with conda..."
conda install -y pysocks=1.7.1=py37_1 || echo "Conda conflict: pysocks=1.7.1=py37_1" >> conda_conflicts.txt
echo "Installing python=3.7.15=h7a1cb2a_1 with conda..."
conda install -y python=3.7.15=h7a1cb2a_1 || echo "Conda conflict: python=3.7.15=h7a1cb2a_1" >> conda_conflicts.txt
echo "Installing python-dateutil=2.8.2=pyhd3eb1b0_0 with conda..."
conda install -y python-dateutil=2.8.2=pyhd3eb1b0_0 || echo "Conda conflict: python-dateutil=2.8.2=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing pytorch=1.13.0=py3.7_cuda11.6_cudnn8.3.2_0 with conda..."
conda install -y pytorch=1.13.0=py3.7_cuda11.6_cudnn8.3.2_0 || echo "Conda conflict: pytorch=1.13.0=py3.7_cuda11.6_cudnn8.3.2_0" >> conda_conflicts.txt
echo "Installing pytorch-cuda=11.6=h867d48c_1 with conda..."
conda install -y pytorch-cuda=11.6=h867d48c_1 || echo "Conda conflict: pytorch-cuda=11.6=h867d48c_1" >> conda_conflicts.txt
echo "Installing pytorch-mutex=1.0=cuda with conda..."
conda install -y pytorch-mutex=1.0=cuda || echo "Conda conflict: pytorch-mutex=1.0=cuda" >> conda_conflicts.txt
echo "Installing pytz=2022.1=py37h06a4308_0 with conda..."
conda install -y pytz=2022.1=py37h06a4308_0 || echo "Conda conflict: pytz=2022.1=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing pyzmq=23.2.0=py37h6a678d5_0 with conda..."
conda install -y pyzmq=23.2.0=py37h6a678d5_0 || echo "Conda conflict: pyzmq=23.2.0=py37h6a678d5_0" >> conda_conflicts.txt
echo "Installing rdkit=2020.09.1.0=py37hd50e099_1 with conda..."
conda install -y rdkit=2020.09.1.0=py37hd50e099_1 || echo "Conda conflict: rdkit=2020.09.1.0=py37hd50e099_1" >> conda_conflicts.txt
echo "Installing readline=8.2=h5eee18b_0 with conda..."
conda install -y readline=8.2=h5eee18b_0 || echo "Conda conflict: readline=8.2=h5eee18b_0" >> conda_conflicts.txt
echo "Installing requests=2.28.1=py37h06a4308_0 with conda..."
conda install -y requests=2.28.1=py37h06a4308_0 || echo "Conda conflict: requests=2.28.1=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing setuptools=65.5.0=py37h06a4308_0 with conda..."
conda install -y setuptools=65.5.0=py37h06a4308_0 || echo "Conda conflict: setuptools=65.5.0=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing six=1.16.0=pyhd3eb1b0_1 with conda..."
conda install -y six=1.16.0=pyhd3eb1b0_1 || echo "Conda conflict: six=1.16.0=pyhd3eb1b0_1" >> conda_conflicts.txt
echo "Installing sqlite=3.40.0=h5082296_0 with conda..."
conda install -y sqlite=3.40.0=h5082296_0 || echo "Conda conflict: sqlite=3.40.0=h5082296_0" >> conda_conflicts.txt
echo "Installing tk=8.6.12=h1ccaba5_0 with conda..."
conda install -y tk=8.6.12=h1ccaba5_0 || echo "Conda conflict: tk=8.6.12=h1ccaba5_0" >> conda_conflicts.txt
echo "Installing torchaudio=0.13.0=py37_cu116 with conda..."
conda install -y torchaudio=0.13.0=py37_cu116 || echo "Conda conflict: torchaudio=0.13.0=py37_cu116" >> conda_conflicts.txt
echo "Installing torchvision=0.14.0=py37_cu116 with conda..."
conda install -y torchvision=0.14.0=py37_cu116 || echo "Conda conflict: torchvision=0.14.0=py37_cu116" >> conda_conflicts.txt
echo "Installing tornado=6.2=py37h5eee18b_0 with conda..."
conda install -y tornado=6.2=py37h5eee18b_0 || echo "Conda conflict: tornado=6.2=py37h5eee18b_0" >> conda_conflicts.txt
echo "Installing traitlets=5.7.1=py37h06a4308_0 with conda..."
conda install -y traitlets=5.7.1=py37h06a4308_0 || echo "Conda conflict: traitlets=5.7.1=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing typing_extensions=4.4.0=py37h06a4308_0 with conda..."
conda install -y typing_extensions=4.4.0=py37h06a4308_0 || echo "Conda conflict: typing_extensions=4.4.0=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing urllib3=1.26.13=py37h06a4308_0 with conda..."
conda install -y urllib3=1.26.13=py37h06a4308_0 || echo "Conda conflict: urllib3=1.26.13=py37h06a4308_0" >> conda_conflicts.txt
echo "Installing wcwidth=0.2.5=pyhd3eb1b0_0 with conda..."
conda install -y wcwidth=0.2.5=pyhd3eb1b0_0 || echo "Conda conflict: wcwidth=0.2.5=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing wheel=0.37.1=pyhd3eb1b0_0 with conda..."
conda install -y wheel=0.37.1=pyhd3eb1b0_0 || echo "Conda conflict: wheel=0.37.1=pyhd3eb1b0_0" >> conda_conflicts.txt
echo "Installing xz=5.2.8=h5eee18b_0 with conda..."
conda install -y xz=5.2.8=h5eee18b_0 || echo "Conda conflict: xz=5.2.8=h5eee18b_0" >> conda_conflicts.txt
echo "Installing zeromq=4.3.4=h2531618_0 with conda..."
conda install -y zeromq=4.3.4=h2531618_0 || echo "Conda conflict: zeromq=4.3.4=h2531618_0" >> conda_conflicts.txt
echo "Installing zlib=1.2.13=h5eee18b_0 with conda..."
conda install -y zlib=1.2.13=h5eee18b_0 || echo "Conda conflict: zlib=1.2.13=h5eee18b_0" >> conda_conflicts.txt
echo "Installing zstd=1.5.2=ha4553b6_0 with conda..."
conda install -y zstd=1.5.2=ha4553b6_0 || echo "Conda conflict: zstd=1.5.2=ha4553b6_0" >> conda_conflicts.txt

# Pip installation
echo "Installing absl-py==1.3.0 with pip..."
pip install absl-py==1.3.0 || echo "Pip conflict: absl-py==1.3.0" >> pip_conflicts.txt
echo "Installing aiohttp==3.8.3 with pip..."
pip install aiohttp==3.8.3 || echo "Pip conflict: aiohttp==3.8.3" >> pip_conflicts.txt
echo "Installing aiosignal==1.3.1 with pip..."
pip install aiosignal==1.3.1 || echo "Pip conflict: aiosignal==1.3.1" >> pip_conflicts.txt
echo "Installing altair==4.2.0 with pip..."
pip install altair==4.2.0 || echo "Pip conflict: altair==4.2.0" >> pip_conflicts.txt
echo "Installing async-timeout==4.0.2 with pip..."
pip install async-timeout==4.0.2 || echo "Pip conflict: async-timeout==4.0.2" >> pip_conflicts.txt
echo "Installing asynctest==0.13.0 with pip..."
pip install asynctest==0.13.0 || echo "Pip conflict: asynctest==0.13.0" >> pip_conflicts.txt
echo "Installing attrs==22.1.0 with pip..."
pip install attrs==22.1.0 || echo "Pip conflict: attrs==22.1.0" >> pip_conflicts.txt
echo "Installing autocommand==2.2.2 with pip..."
pip install autocommand==2.2.2 || echo "Pip conflict: autocommand==2.2.2" >> pip_conflicts.txt
echo "Installing backports-zoneinfo==0.2.1 with pip..."
pip install backports-zoneinfo==0.2.1 || echo "Pip conflict: backports-zoneinfo==0.2.1" >> pip_conflicts.txt
echo "Installing beautifulsoup4==4.11.1 with pip..."
pip install beautifulsoup4==4.11.1 || echo "Pip conflict: beautifulsoup4==4.11.1" >> pip_conflicts.txt
echo "Installing bertviz==1.4.0 with pip..."
pip install bertviz==1.4.0 || echo "Pip conflict: bertviz==1.4.0" >> pip_conflicts.txt
echo "Installing blinker==1.5 with pip..."
pip install blinker==1.5 || echo "Pip conflict: blinker==1.5" >> pip_conflicts.txt
echo "Installing boto3==1.26.29 with pip..."
pip install boto3==1.26.29 || echo "Pip conflict: boto3==1.26.29" >> pip_conflicts.txt
echo "Installing botocore==1.29.29 with pip..."
pip install botocore==1.29.29 || echo "Pip conflict: botocore==1.29.29" >> pip_conflicts.txt
echo "Installing cachetools==5.2.0 with pip..."
pip install cachetools==5.2.0 || echo "Pip conflict: cachetools==5.2.0" >> pip_conflicts.txt
echo "Installing cheroot==9.0.0 with pip..."
pip install cheroot==9.0.0 || echo "Pip conflict: cheroot==9.0.0" >> pip_conflicts.txt
echo "Installing cherrypy==18.8.0 with pip..."
pip install cherrypy==18.8.0 || echo "Pip conflict: cherrypy==18.8.0" >> pip_conflicts.txt
echo "Installing click==8.1.3 with pip..."
pip install click==8.1.3 || echo "Pip conflict: click==8.1.3" >> pip_conflicts.txt
echo "Installing cloudpickle==2.2.0 with pip..."
pip install cloudpickle==2.2.0 || echo "Pip conflict: cloudpickle==2.2.0" >> pip_conflicts.txt
echo "Installing colour==0.1.5 with pip..."
pip install colour==0.1.5 || echo "Pip conflict: colour==0.1.5" >> pip_conflicts.txt
echo "Installing commonmark==0.9.1 with pip..."
pip install commonmark==0.9.1 || echo "Pip conflict: commonmark==0.9.1" >> pip_conflicts.txt
echo "Installing cycler==0.11.0 with pip..."
pip install cycler==0.11.0 || echo "Pip conflict: cycler==0.11.0" >> pip_conflicts.txt
echo "Installing datasets==2.7.1 with pip..."
pip install datasets==2.7.1 || echo "Pip conflict: datasets==2.7.1" >> pip_conflicts.txt
echo "Installing dill==0.3.6 with pip..."
pip install dill==0.3.6 || echo "Pip conflict: dill==0.3.6" >> pip_conflicts.txt
echo "Installing docker-pycreds==0.4.0 with pip..."
pip install docker-pycreds==0.4.0 || echo "Pip conflict: docker-pycreds==0.4.0" >> pip_conflicts.txt
echo "Installing dominate==2.7.0 with pip..."
pip install dominate==2.7.0 || echo "Pip conflict: dominate==2.7.0" >> pip_conflicts.txt
echo "Installing drfp==0.3.2 with pip..."
pip install drfp==0.3.2 || echo "Pip conflict: drfp==0.3.2" >> pip_conflicts.txt
echo "Installing et-xmlfile==1.1.0 with pip..."
pip install et-xmlfile==1.1.0 || echo "Pip conflict: et-xmlfile==1.1.0" >> pip_conflicts.txt
echo "Installing faerun==0.3.20 with pip..."
pip install faerun==0.3.20 || echo "Pip conflict: faerun==0.3.20" >> pip_conflicts.txt
echo "Installing filelock==3.8.2 with pip..."
pip install filelock==3.8.2 || echo "Pip conflict: filelock==3.8.2" >> pip_conflicts.txt
echo "Installing flask==2.0.2 with pip..."
pip install flask==2.0.2 || echo "Pip conflict: flask==2.0.2" >> pip_conflicts.txt
echo "Installing flask-bootstrap==3.3.7.1 with pip..."
pip install flask-bootstrap==3.3.7.1 || echo "Pip conflict: flask-bootstrap==3.3.7.1" >> pip_conflicts.txt
echo "Installing flask-wtf==1.0.0 with pip..."
pip install flask-wtf==1.0.0 || echo "Pip conflict: flask-wtf==1.0.0" >> pip_conflicts.txt
echo "Installing frozenlist==1.3.3 with pip..."
pip install frozenlist==1.3.3 || echo "Pip conflict: frozenlist==1.3.3" >> pip_conflicts.txt
echo "Installing fsspec==2022.11.0 with pip..."
pip install fsspec==2022.11.0 || echo "Pip conflict: fsspec==2022.11.0" >> pip_conflicts.txt
echo "Installing future==0.18.3 with pip..."
pip install future==0.18.3 || echo "Pip conflict: future==0.18.3" >> pip_conflicts.txt
echo "Installing gdown==4.6.0 with pip..."
pip install gdown==4.6.0 || echo "Pip conflict: gdown==4.6.0" >> pip_conflicts.txt
echo "Installing gitdb==4.0.10 with pip..."
pip install gitdb==4.0.10 || echo "Pip conflict: gitdb==4.0.10" >> pip_conflicts.txt
echo "Installing gitpython==3.1.29 with pip..."
pip install gitpython==3.1.29 || echo "Pip conflict: gitpython==3.1.29" >> pip_conflicts.txt
echo "Installing google-auth==2.15.0 with pip..."
pip install google-auth==2.15.0 || echo "Pip conflict: google-auth==2.15.0" >> pip_conflicts.txt
echo "Installing google-auth-oauthlib==0.4.6 with pip..."
pip install google-auth-oauthlib==0.4.6 || echo "Pip conflict: google-auth-oauthlib==0.4.6" >> pip_conflicts.txt
echo "Installing grpcio==1.51.1 with pip..."
pip install grpcio==1.51.1 || echo "Pip conflict: grpcio==1.51.1" >> pip_conflicts.txt
echo "Installing huggingface-hub==0.11.1 with pip..."
pip install huggingface-hub==0.11.1 || echo "Pip conflict: huggingface-hub==0.11.1" >> pip_conflicts.txt
echo "Installing hyperopt==0.2.7 with pip..."
pip install hyperopt==0.2.7 || echo "Pip conflict: hyperopt==0.2.7" >> pip_conflicts.txt
echo "Installing importlib-metadata==5.1.0 with pip..."
pip install importlib-metadata==5.1.0 || echo "Pip conflict: importlib-metadata==5.1.0" >> pip_conflicts.txt
echo "Installing importlib-resources==5.10.1 with pip..."
pip install importlib-resources==5.10.1 || echo "Pip conflict: importlib-resources==5.10.1" >> pip_conflicts.txt
echo "Installing inflect==6.0.2 with pip..."
pip install inflect==6.0.2 || echo "Pip conflict: inflect==6.0.2" >> pip_conflicts.txt
echo "Installing ipython==7.34.0 with pip..."
pip install ipython==7.34.0 || echo "Pip conflict: ipython==7.34.0" >> pip_conflicts.txt
echo "Installing itsdangerous==2.1.2 with pip..."
pip install itsdangerous==2.1.2 || echo "Pip conflict: itsdangerous==2.1.2" >> pip_conflicts.txt
echo "Installing jaraco-classes==3.2.3 with pip..."
pip install jaraco-classes==3.2.3 || echo "Pip conflict: jaraco-classes==3.2.3" >> pip_conflicts.txt
echo "Installing jaraco-collections==3.8.0 with pip..."
pip install jaraco-collections==3.8.0 || echo "Pip conflict: jaraco-collections==3.8.0" >> pip_conflicts.txt
echo "Installing jaraco-context==4.2.0 with pip..."
pip install jaraco-context==4.2.0 || echo "Pip conflict: jaraco-context==4.2.0" >> pip_conflicts.txt
echo "Installing jaraco-functools==3.5.2 with pip..."
pip install jaraco-functools==3.5.2 || echo "Pip conflict: jaraco-functools==3.5.2" >> pip_conflicts.txt
echo "Installing jaraco-text==3.11.0 with pip..."
pip install jaraco-text==3.11.0 || echo "Pip conflict: jaraco-text==3.11.0" >> pip_conflicts.txt
echo "Installing jedi==0.18.2 with pip..."
pip install jedi==0.18.2 || echo "Pip conflict: jedi==0.18.2" >> pip_conflicts.txt
echo "Installing jinja2==3.1.2 with pip..."
pip install jinja2==3.1.2 || echo "Pip conflict: jinja2==3.1.2" >> pip_conflicts.txt
echo "Installing jmespath==1.0.1 with pip..."
pip install jmespath==1.0.1 || echo "Pip conflict: jmespath==1.0.1" >> pip_conflicts.txt
echo "Installing joblib==1.2.0 with pip..."
pip install joblib==1.2.0 || echo "Pip conflict: joblib==1.2.0" >> pip_conflicts.txt
echo "Installing jsonschema==4.17.3 with pip..."
pip install jsonschema==4.17.3 || echo "Pip conflict: jsonschema==4.17.3" >> pip_conflicts.txt
echo "Installing kiwisolver==1.4.4 with pip..."
pip install kiwisolver==1.4.4 || echo "Pip conflict: kiwisolver==1.4.4" >> pip_conflicts.txt
echo "Installing llvmlite==0.39.1 with pip..."
pip install llvmlite==0.39.1 || echo "Pip conflict: llvmlite==0.39.1" >> pip_conflicts.txt
echo "Installing markdown==3.4.1 with pip..."
pip install markdown==3.4.1 || echo "Pip conflict: markdown==3.4.1" >> pip_conflicts.txt
echo "Installing markupsafe==2.1.1 with pip..."
pip install markupsafe==2.1.1 || echo "Pip conflict: markupsafe==2.1.1" >> pip_conflicts.txt
echo "Installing matplotlib==3.2.2 with pip..."
pip install matplotlib==3.2.2 || echo "Pip conflict: matplotlib==3.2.2" >> pip_conflicts.txt
echo "Installing more-itertools==9.0.0 with pip..."
pip install more-itertools==9.0.0 || echo "Pip conflict: more-itertools==9.0.0" >> pip_conflicts.txt
echo "Installing multidict==6.0.3 with pip..."
pip install multidict==6.0.3 || echo "Pip conflict: multidict==6.0.3" >> pip_conflicts.txt
echo "Installing multiprocess==0.70.14 with pip..."
pip install multiprocess==0.70.14 || echo "Pip conflict: multiprocess==0.70.14" >> pip_conflicts.txt
echo "Installing networkx==2.6.3 with pip..."
pip install networkx==2.6.3 || echo "Pip conflict: networkx==2.6.3" >> pip_conflicts.txt
echo "Installing numba==0.56.4 with pip..."
pip install numba==0.56.4 || echo "Pip conflict: numba==0.56.4" >> pip_conflicts.txt
echo "Installing oauthlib==3.2.2 with pip..."
pip install oauthlib==3.2.2 || echo "Pip conflict: oauthlib==3.2.2" >> pip_conflicts.txt
echo "Installing openpyxl==3.1.0 with pip..."
pip install openpyxl==3.1.0 || echo "Pip conflict: openpyxl==3.1.0" >> pip_conflicts.txt
echo "Installing pandarallel==1.6.3 with pip..."
pip install pandarallel==1.6.3 || echo "Pip conflict: pandarallel==1.6.3" >> pip_conflicts.txt
echo "Installing pathtools==0.1.2 with pip..."
pip install pathtools==0.1.2 || echo "Pip conflict: pathtools==0.1.2" >> pip_conflicts.txt
echo "Installing pkgutil-resolve-name==1.3.10 with pip..."
pip install pkgutil-resolve-name==1.3.10 || echo "Pip conflict: pkgutil-resolve-name==1.3.10" >> pip_conflicts.txt
echo "Installing portend==3.1.0 with pip..."
pip install portend==3.1.0 || echo "Pip conflict: portend==3.1.0" >> pip_conflicts.txt
echo "Installing promise==2.3 with pip..."
pip install promise==2.3 || echo "Pip conflict: promise==2.3" >> pip_conflicts.txt
echo "Installing prompt-toolkit==3.0.36 with pip..."
pip install prompt-toolkit==3.0.36 || echo "Pip conflict: prompt-toolkit==3.0.36" >> pip_conflicts.txt
echo "Installing protobuf==3.20.1 with pip..."
pip install protobuf==3.20.1 || echo "Pip conflict: protobuf==3.20.1" >> pip_conflicts.txt
echo "Installing psutil==5.9.4 with pip..."
pip install psutil==5.9.4 || echo "Pip conflict: psutil==5.9.4" >> pip_conflicts.txt
echo "Installing py4j==0.10.9.7 with pip..."
pip install py4j==0.10.9.7 || echo "Pip conflict: py4j==0.10.9.7" >> pip_conflicts.txt
echo "Installing pyarrow==10.0.1 with pip..."
pip install pyarrow==10.0.1 || echo "Pip conflict: pyarrow==10.0.1" >> pip_conflicts.txt
echo "Installing pyasn1==0.4.8 with pip..."
pip install pyasn1==0.4.8 || echo "Pip conflict: pyasn1==0.4.8" >> pip_conflicts.txt
echo "Installing pyasn1-modules==0.2.8 with pip..."
pip install pyasn1-modules==0.2.8 || echo "Pip conflict: pyasn1-modules==0.2.8" >> pip_conflicts.txt
echo "Installing pydantic==1.10.2 with pip..."
pip install pydantic==1.10.2 || echo "Pip conflict: pydantic==1.10.2" >> pip_conflicts.txt
echo "Installing pydeck==0.8.0 with pip..."
pip install pydeck==0.8.0 || echo "Pip conflict: pydeck==0.8.0" >> pip_conflicts.txt
echo "Installing pygments==2.13.0 with pip..."
pip install pygments==2.13.0 || echo "Pip conflict: pygments==2.13.0" >> pip_conflicts.txt
echo "Installing pympler==1.0.1 with pip..."
pip install pympler==1.0.1 || echo "Pip conflict: pympler==1.0.1" >> pip_conflicts.txt
echo "Installing pyrsistent==0.19.2 with pip..."
pip install pyrsistent==0.19.2 || echo "Pip conflict: pyrsistent==0.19.2" >> pip_conflicts.txt
echo "Installing pytz-deprecation-shim==0.1.0.post0 with pip..."
pip install pytz-deprecation-shim==0.1.0.post0 || echo "Pip conflict: pytz-deprecation-shim==0.1.0.post0" >> pip_conflicts.txt
echo "Installing pyyaml==6.0 with pip..."
pip install pyyaml==6.0 || echo "Pip conflict: pyyaml==6.0" >> pip_conflicts.txt
echo "Installing rdkit-pypi==2022.9.3 with pip..."
pip install rdkit-pypi==2022.9.3 || echo "Pip conflict: rdkit-pypi==2022.9.3" >> pip_conflicts.txt
echo "Installing regex==2022.10.31 with pip..."
pip install regex==2022.10.31 || echo "Pip conflict: regex==2022.10.31" >> pip_conflicts.txt
echo "Installing requests-oauthlib==1.3.1 with pip..."
pip install requests-oauthlib==1.3.1 || echo "Pip conflict: requests-oauthlib==1.3.1" >> pip_conflicts.txt
echo "Installing responses==0.18.0 with pip..."
pip install responses==0.18.0 || echo "Pip conflict: responses==0.18.0" >> pip_conflicts.txt
echo "Installing rich==12.6.0 with pip..."
pip install rich==12.6.0 || echo "Pip conflict: rich==12.6.0" >> pip_conflicts.txt
echo "Installing rsa==4.9 with pip..."
pip install rsa==4.9 || echo "Pip conflict: rsa==4.9" >> pip_conflicts.txt
echo "Installing rxnfp==0.1.0 with pip..."
pip install rxnfp==0.1.0 || echo "Pip conflict: rxnfp==0.1.0" >> pip_conflicts.txt
echo "Installing s3transfer==0.6.0 with pip..."
pip install s3transfer==0.6.0 || echo "Pip conflict: s3transfer==0.6.0" >> pip_conflicts.txt
echo "Installing sacremoses==0.0.53 with pip..."
pip install sacremoses==0.0.53 || echo "Pip conflict: sacremoses==0.0.53" >> pip_conflicts.txt
echo "Installing scikit-learn==0.23.1 with pip..."
pip install scikit-learn==0.23.1 || echo "Pip conflict: scikit-learn==0.23.1" >> pip_conflicts.txt
echo "Installing scipy==1.4.1 with pip..."
pip install scipy==1.4.1 || echo "Pip conflict: scipy==1.4.1" >> pip_conflicts.txt
echo "Installing seaborn==0.12.1 with pip..."
pip install seaborn==0.12.1 || echo "Pip conflict: seaborn==0.12.1" >> pip_conflicts.txt
echo "Installing semver==2.13.0 with pip..."
pip install semver==2.13.0 || echo "Pip conflict: semver==2.13.0" >> pip_conflicts.txt
echo "Installing sentencepiece==0.1.97 with pip..."
pip install sentencepiece==0.1.97 || echo "Pip conflict: sentencepiece==0.1.97" >> pip_conflicts.txt
echo "Installing sentry-sdk==1.11.1 with pip..."
pip install sentry-sdk==1.11.1 || echo "Pip conflict: sentry-sdk==1.11.1" >> pip_conflicts.txt
echo "Installing seqeval==1.2.2 with pip..."
pip install seqeval==1.2.2 || echo "Pip conflict: seqeval==1.2.2" >> pip_conflicts.txt
echo "Installing setproctitle==1.3.2 with pip..."
pip install setproctitle==1.3.2 || echo "Pip conflict: setproctitle==1.3.2" >> pip_conflicts.txt
echo "Installing shap==0.41.0 with pip..."
pip install shap==0.41.0 || echo "Pip conflict: shap==0.41.0" >> pip_conflicts.txt
echo "Installing shortuuid==1.0.11 with pip..."
pip install shortuuid==1.0.11 || echo "Pip conflict: shortuuid==1.0.11" >> pip_conflicts.txt
echo "Installing simpletransformers==0.61.13 with pip..."
pip install simpletransformers==0.61.13 || echo "Pip conflict: simpletransformers==0.61.13" >> pip_conflicts.txt
echo "Installing slicer==0.0.7 with pip..."
pip install slicer==0.0.7 || echo "Pip conflict: slicer==0.0.7" >> pip_conflicts.txt
echo "Installing smmap==5.0.0 with pip..."
pip install smmap==5.0.0 || echo "Pip conflict: smmap==5.0.0" >> pip_conflicts.txt
echo "Installing soupsieve==2.3.2.post1 with pip..."
pip install soupsieve==2.3.2.post1 || echo "Pip conflict: soupsieve==2.3.2.post1" >> pip_conflicts.txt
echo "Installing streamlit==1.15.2 with pip..."
pip install streamlit==1.15.2 || echo "Pip conflict: streamlit==1.15.2" >> pip_conflicts.txt
echo "Installing tempora==5.1.0 with pip..."
pip install tempora==5.1.0 || echo "Pip conflict: tempora==5.1.0" >> pip_conflicts.txt
echo "Installing tensorboard==2.11.0 with pip..."
pip install tensorboard==2.11.0 || echo "Pip conflict: tensorboard==2.11.0" >> pip_conflicts.txt
echo "Installing tensorboard-data-server==0.6.1 with pip..."
pip install tensorboard-data-server==0.6.1 || echo "Pip conflict: tensorboard-data-server==0.6.1" >> pip_conflicts.txt
echo "Installing tensorboard-plugin-wit==1.8.1 with pip..."
pip install tensorboard-plugin-wit==1.8.1 || echo "Pip conflict: tensorboard-plugin-wit==1.8.1" >> pip_conflicts.txt
echo "Installing tensorboardx==2.5.1 with pip..."
pip install tensorboardx==2.5.1 || echo "Pip conflict: tensorboardx==2.5.1" >> pip_conflicts.txt
echo "Installing threadpoolctl==3.1.0 with pip..."
pip install threadpoolctl==3.1.0 || echo "Pip conflict: threadpoolctl==3.1.0" >> pip_conflicts.txt
echo "Installing tokenizers==0.12.1 with pip..."
pip install tokenizers==0.12.1 || echo "Pip conflict: tokenizers==0.12.1" >> pip_conflicts.txt
echo "Installing toml==0.10.2 with pip..."
pip install toml==0.10.2 || echo "Pip conflict: toml==0.10.2" >> pip_conflicts.txt
echo "Installing toolz==0.12.0 with pip..."
pip install toolz==0.12.0 || echo "Pip conflict: toolz==0.12.0" >> pip_conflicts.txt
echo "Installing tqdm==4.64.1 with pip..."
pip install tqdm==4.64.1 || echo "Pip conflict: tqdm==4.64.1" >> pip_conflicts.txt
echo "Installing transformers==4.18.0 with pip..."
pip install transformers==4.18.0 || echo "Pip conflict: transformers==4.18.0" >> pip_conflicts.txt
echo "Installing tzdata==2022.7 with pip..."
pip install tzdata==2022.7 || echo "Pip conflict: tzdata==2022.7" >> pip_conflicts.txt
echo "Installing tzlocal==4.2 with pip..."
pip install tzlocal==4.2 || echo "Pip conflict: tzlocal==4.2" >> pip_conflicts.txt
echo "Installing ujson==5.6.0 with pip..."
pip install ujson==5.6.0 || echo "Pip conflict: ujson==5.6.0" >> pip_conflicts.txt
echo "Installing validators==0.20.0 with pip..."
pip install validators==0.20.0 || echo "Pip conflict: validators==0.20.0" >> pip_conflicts.txt
echo "Installing visitor==0.1.3 with pip..."
pip install visitor==0.1.3 || echo "Pip conflict: visitor==0.1.3" >> pip_conflicts.txt
echo "Installing wandb==0.13.6 with pip..."
pip install wandb==0.13.6 || echo "Pip conflict: wandb==0.13.6" >> pip_conflicts.txt
echo "Installing watchdog==2.2.0 with pip..."
pip install watchdog==2.2.0 || echo "Pip conflict: watchdog==2.2.0" >> pip_conflicts.txt
echo "Installing werkzeug==2.2.2 with pip..."
pip install werkzeug==2.2.2 || echo "Pip conflict: werkzeug==2.2.2" >> pip_conflicts.txt
echo "Installing wtforms==3.0.1 with pip..."
pip install wtforms==3.0.1 || echo "Pip conflict: wtforms==3.0.1" >> pip_conflicts.txt
echo "Installing xgboost==1.6.2 with pip..."
pip install xgboost==1.6.2 || echo "Pip conflict: xgboost==1.6.2" >> pip_conflicts.txt
echo "Installing xlwt==1.3.0 with pip..."
pip install xlwt==1.3.0 || echo "Pip conflict: xlwt==1.3.0" >> pip_conflicts.txt
echo "Installing xxhash==3.1.0 with pip..."
pip install xxhash==3.1.0 || echo "Pip conflict: xxhash==3.1.0" >> pip_conflicts.txt
echo "Installing yarl==1.8.2 with pip..."
pip install yarl==1.8.2 || echo "Pip conflict: yarl==1.8.2" >> pip_conflicts.txt
echo "Installing zc-lockfile==2.0 with pip..."
pip install zc-lockfile==2.0 || echo "Pip conflict: zc-lockfile==2.0" >> pip_conflicts.txt
echo "Installing zipp==3.11.0 with pip..."
pip install zipp==3.11.0 || echo "Pip conflict: zipp==3.11.0" >> pip_conflicts.txt
