U=2 Example - Non-PH symmetric
==============================
This is an example that is not Particle-Hole symmetric. These files illustrate the subtle input differences between the other examples and non-symmetric cases. For full details about the output of `Maxent`, please see the other PDFs in the example folder. 

### Files
**_G_im_** = Im[G] 
  - Column Format: `iw_n Im[G]_n error_n  `

**_G_re_** = Re[G] 
  - Column Format: `iw_n Re[G]_n error_n` 

**_Gomegain_** = Input format of G for `Maxent`
  - Column Format: `iwn_n Re[G]_n error_re_n Im[G]_n error_im_n`  

**_G_tau_** = G(tau)
  - Column Format `tau_n G(tau)_n error_n`  

**_Selfenergy_** = Sigma(iw_n)  
  - Column Format `iwn_n Re[Sigma] error__rn Im[Sigma] error_i_n

**_frequency.param_** = parameter file for Matsubara frequency data  
**_time.param_** = parameter file for time data  
**_self.param_** = parameter file for self-energy data
