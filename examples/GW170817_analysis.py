"""
This script will analyse the open data for GW170817, as well as the associated
GRB and optical data to jointly determine the parameters of all the models.
"""

from __future__ import division
import bilby
import multimediabns

outdir = 'outdir'
label = 'GW170817'
time_of_event = bilby.gw.utils.get_event_time(label)

# set up the logger for the output from BILBY
bilby.core.utils.setup_logger(outdir=outdir, label=label)

# ====================================================================
# ========= CURRENTLY CANNOT USE LOSC TO GET THE EVENT DATA & PSD ====
# ====================================================================
# interferometers = bilby.gw.detector.get_event_data(label, tag='CLN',
#     outdir=outdir, psd_offset=-1024,
#     psd_duration=100, plot=True)


# ====================================================================
# ========= FOR NOW, SET UP INJECTION INTO aLIGO =====================
# ====================================================================
GW_params = dict(
    luminosity_distance=80.,
    chirp_mass=1.1188, mass_ratio=0.9, chi_1=0.00, chi_2=0.00,
    iota=0.4, psi=2.659, phase=1.3, geocent_time=time_of_event,
    ra=1.375, dec=-1.2108, lambda_1=400, lambda_2=400)

duration = 8.
sampling_frequency = 2 * 1570.
start_time = GW_params['geocent_time'] + 2 - duration
minimum_frequency = 40.0

# Set up the interferometers
interferometers = bilby.gw.detector.InterferometerList(['H1', 'L1', 'V1'])
for interferometer in interferometers:
    interferometer.minimum_frequency = minimum_frequency

interferometers.set_strain_data_from_power_spectral_densities(
        sampling_frequency=sampling_frequency, duration=duration,
        start_time=start_time)

# WAVEFORM GENERATOR
waveform_arguments = dict(waveform_approximant='TaylorF2',
                          reference_frequency=50.0, minimum_frequency=minimum_frequency)
# Create the waveform_generator using a LAL Binary Neutron Star source function
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_neutron_star,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
    waveform_arguments=waveform_arguments)

# inject the signal with GW_params into the detector network
interferometers.inject_signal(parameters=GW_params,
                              waveform_generator=waveform_generator)

# PRIORS
GW_priors = bilby.gw.prior.BNSPriorDict()
GW_priors.pop('mass_1')
GW_priors.pop('mass_2')
GW_priors.pop('lambda_1')
GW_priors.pop('lambda_2')
GW_priors['chirp_mass'] = bilby.core.prior.Uniform(
    0.87, 1.74, name='chirp_mass', unit='$M_{\\odot}$')
GW_priors['mass_ratio'] = bilby.core.prior.Uniform(
    0.5, 1.0, name='mass_ratio')
GW_priors['lambda_tilde'] = bilby.core.prior.Uniform(0, 5000, name='lambda_tilde')
GW_priors['delta_lambda'] = bilby.core.prior.Uniform(
    -5000, 5000, name='delta_lambda')

GW_priors['geocent_time'] = bilby.core.prior.Uniform(
    minimum=GW_params['geocent_time'] - 1,
    maximum=GW_params['geocent_time'] + 1,
    name='geocent_time', latex_label='$t_c$', unit='$s$')

for key in ['chi_1', 'chi_2']:
    GW_priors[key] = GW_params[key]

# LIKELIHOOD:
GW_likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=interferometers, waveform_generator=waveform_generator,
    time_marginalization=False, phase_marginalization=False,
    distance_marginalization=False, priors=GW_priors)

# Note that for when we want to add in the multiple likelihoods, we simply
# need the line:
# joint_likelihood = bilby.core.JointLikielihood(GW_likelihood, GRB_likelihood...)
# would also need to remove the conversion function in the result section of code

# Run sampler
result = bilby.run_sampler(
    likelihood=GW_likelihood, priors=GW_priors, sampler='dynesty', npoints=n_live,
    injection_parameters=GW_params, outdir=outdir, label=label,
    conversion_function=bilby.gw.conversion.generate_all_bns_parameters,
    walks=60)

result.plot_corner()
