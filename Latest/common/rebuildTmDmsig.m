function [sig_lms, sig_rls, sig_lms_dd, sig_rls_dd] = rebuildTmDmsig(bitsPerSym, sps, eBW, y_lms, y_rls, y_lms_dd, y_rls_dd)
rClean_lms = (real(y_lms) >0);
rClean_rls = (real(y_rls) >0);
rClean_lms_dd = (real(y_lms_dd) >0);
rClean_rls_dd = (real(y_rls_dd) >0);
sig_lms = psk_mod(bitsPerSym, sps, eBW, rClean_lms);
sig_rls= psk_mod(bitsPerSym, sps, eBW, rClean_rls);
sig_lms_dd = psk_mod(bitsPerSym, sps, eBW, rClean_lms_dd);
sig_rls_dd = psk_mod(bitsPerSym, sps, eBW, rClean_rls_dd);
end

