function smb = calc_smb_gradient(hSn, zela, grad_smb)

smb = grad_smb*(hSn - zela);

max_smb = 0.5;
min_smb = -3.0;

smb(smb>max_smb) = max_smb;
smb(smb<min_smb) = min_smb;

end