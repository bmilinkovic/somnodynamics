function [dopt,Lopt,conv,iopt,sopt,cput,ohist] = canonical_opt_ddx_mruns( ...
    CAK, L0, niters, gdes, gdsig0, gdls, gdtol, hist, pp, stage_label)
% Canonical pipeline wrapper for proxy-DDX multi-restart optimisation.
%
% This preserves the original two-stage SSDI structure:
%   pre-stage  -> ddx proxy objective (CAK-based)
%   final-stage -> canonical spectral objective (Eq. 25/26; see canonical_opt_dds_mruns)

if nargin < 10 || isempty(stage_label), stage_label = 'pre'; end %#ok<NASGU>

[dopt,Lopt,conv,iopt,sopt,cput,ohist] = opt_gd_ddx_mruns( ...
    CAK, L0, niters, gdes, gdsig0, gdls, gdtol, hist, pp);

end
