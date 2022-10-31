# Estimation code for "Incorporating Time-Varying Covariates in a Simple Mixture Model for Discrete-Time Duration Data"

Pete Fader and Bruce Hardie provide a technical note <https://brucehardie.com/notes/037/> on incorporating time-varying covariates for a discrete-time duration mixture model.  It is an extention to the (shifted) Beta-Geometric Distribution, now with covariates.  The authors call it G2G+covariates model.

However, a code on how to estimate this model is not available.  

I created the MLE estimation code using R with simulated data.  

G2G_varying_MLE.R is the version that matches the technical note.  

G2G_static_MLE.R is a simpler version with static covariates.  I first created it for practice.  It is uploaded to this respostory nevertheless.  If you are trying to learn, it may be easier to go thru this code first.

Please contact me (kaloklee@gmail.com) if you have any questions or comments.  

Some references for the (shifted) Beta-Geometric Distribution:

Fader, Peter S. and Bruce G.S. Hardie (2007), "How to Project Customer Retention," Journal of Interactive Marketing, 31 (1), 76-90.

Lee, Ka Lok, Peter S. Fader, and Bruce G.S. Hardie (2007), “How to Project Patient Persistency,” Foresight: The International Journal of Applied Forecasting, Issue 8, 31-35. 

