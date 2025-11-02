function smb = calc_smb_gradient(hSn, zela, grad_smb)

smb = grad_smb*(hSn - zela);

max_smb = 2.0;
min_smb = -0.2;

smb(smb>max_smb) = max_smb;
smb(smb<min_smb) = min_smb;

end