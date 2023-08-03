# RFI-suppressed_L-band_SMOS_multi-angular_TB_Refinement

Source codes that were used to produce the dataset associated with the following research:

Title: An RFI-suppressed SMOS L-band multi-angular brightness temperature dataset spanning over a decade (since 2010) Authors: Zhiqing Peng, Tianjie Zhao, Jiancheng Shi, Yann H. Kerr, Nemesio J. Rodríguez-Fernández, Panpan Yao, Tao Che.

The RFI-suppressed SMOS multi-angular brightness temperature can be obtained from the SMOS SCLF1C product (version 724) by “two_step_refinement_main.exe” (unzip from the 5 rar files). 

For example, after unzipping the SCLF1C product file, the ‘.HDR’ and the ‘.DBL’ format file should be in the same folder, then users can use it with the following command:

two_step_refinement_main.exe –infile SM_REPR_MIR_SCLF1C_20160601T020556_20160601T023225_724_200_1.DBL(DBL data full path)

The output file format of “two_step_refinement_main.exe” is hdf5 and in the projection of ISEA4H9 15-KM, users can convert it to netCDF4 with the script named “convert_twostep_h5f2nc.py”.

“read_plot_refined_smos_multi-angular_TB.py” contains codes for reading and plotting the refined SMOS multi-angular TB in the projection of ISEA4H9 and EASE-GRID 2.0.

Related dataset:
The associated dataset can be downloaded from TPDC：https://doi.org/10.11888/Terre.tpdc.300406 or https://cstr.cn/18406.11.Terre.tpdc.300406
