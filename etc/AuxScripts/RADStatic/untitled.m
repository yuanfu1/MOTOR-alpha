close all;
clear all;

vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9763_rwnd');
data = vardata(vardata~=0);
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9754_rwnd');
data = [data; vardata(vardata~=0)];
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9751_rwnd');
data = [data; vardata(vardata~=0)];
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9660_rwnd');
data = [data; vardata(vardata~=0)];
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9758_rwnd');
data = [data; vardata(vardata~=0)];
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9762_rwnd');
data = [data; vardata(vardata~=0)];
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9759_rwnd');
data = [data; vardata(vardata~=0)];
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9200_rwnd');
data = [data; vardata(vardata~=0)];
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9755_rwnd');
data = [data; vardata(vardata~=0)];
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9662_rwnd');
data = [data; vardata(vardata~=0)];
vardata = ncread('/Volumes/products/MOTOR/SurfaceAnalysis.alpha/220527_0000/output/MOTOR-3DVar_innovation_G06.nc','RAD_Z9753_rwnd');
data = [data; vardata(vardata~=0)];


[a,b] = hist(data,1000);
std(data) % = 2.10
mean(data) % = -0.03
bar(b, (a))

