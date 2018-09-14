import logging
import numpy as np
from math import sqrt

from ..datamodels import dqflags

import warnings
warnings.simplefilter('ignore', np.RankWarning)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class PixelRamp():
    """
    Base class for a PixelRamp object
    """

    def __init__(self, counts, errs, dqs, times):

        # Copy the inputs to class attributes
        self.counts = counts
        self.errs = errs
        self.dqs = dqs
        self.times = times

        self.group_time = times[1] - times[0]
        self.num_groups = len(counts)

        self.groups_with_CR = []
        self.nodes = np.zeros((1, 2), dtype=int)
        self.nodes[0][1] = self.num_groups
        self.num_semiramps = len(self.nodes)

        # Setup semi-ramps based on any existing CR and SAT flags;
        # if there aren't any existing CR flags, there will be
        # one continuous semi-ramp for the whole range of data
        for group in range(self.num_groups):

            # If we hit a saturation flag, mark the previous group as the
            # end of the valid range for this semi-ramp and break
            if self.dqs[group] & dqflags.group['SATURATED']:
                self.nodes[self.num_semiramps - 1][1] = group
                self.num_groups = group
                break

            # If we hit a CR flag, mark the previous group as the end of the
            # valid range for this semi-ramp and start a new semi-ramp at
            # this group
            if self.dqs[group] & dqflags.group['JUMP_DET']:
                self.groups_with_CR.append(group)
                self.nodes[self.num_semiramps - 1][1] = group
                self.nodes = np.append(self.nodes,
                                       [[group, self.num_groups]], axis=0)
                self.num_semiramps += 1


def find_crs(data, err, gdq, times, read_noise, rejection_threshold,
    signal_threshold, median_slopes):

    # Get the attributes of the input data array
    (nints, ngroups, nrows, ncols) = data.shape

    # Loop over multiple integrations
    for integration in range(nints):

        # Loop over all pixels of this integration
        for row in range(nrows):
            for col in range(ncols):

                # Create a PixelRamp object for this pixel
                counts = data[integration, :, row, col]
                errs = err[integration, :, row, col]
                gdqs = gdq[integration, :, row, col]
                ramp = PixelRamp(counts, errs, gdqs, times)

                # Compute the slope for this pixel
                #slope = mean_ramp_slope (ramp, read_noise[row,col])
                slope = median_slopes[integration, row, col]

                # If this pixel is in the read noise regime,
                # apply the y-intercept method for finding CRs
                if slope < signal_threshold:

                    # Find CRs
                    yint(ramp, read_noise[row, col], rejection_threshold)

                    # Add any new CR flags to input DQ array
                    gdq[integration, :, row, col] = np.bitwise_or(
                        gdq[integration, :, row, col], ramp.dqs)

        # next pixel
    # next integration

    return


def mean_ramp_slope(ramp, read_noise):
    """
    Compute the slope of a entire pixel ramp, using the weighted mean
    slope of all defined semi-ramps.
    """

    sum_slopes = sum_weights = 0.

    # Loop over all semi-ramps, computing slope for each
    for semiramp in range(ramp.num_semiramps):

        # Get the start and end positions for this semi-ramp
        beg = ramp.nodes[semiramp][0]
        end = ramp.nodes[semiramp][1]

        # Check to make sure there are enough points for a fit
        if (end - beg) < 2:
            slope = weight = 0.

        else:

            # Compute the slope and err for this semi-ramp
            coeffs, resid, rank, sing, rcond = \
                np.polyfit(ramp.times[beg:end], ramp.counts[beg:end], 1, full=True)
            # Compute the weight for this slope
            if len(resid) == 0:
                slope = weight = 0.
            else:
                slope = coeffs[0]
                weight = 1. / resid

        # Accumulate the weighted slope
        sum_slopes += slope * weight
        sum_weights += weight

    # Return the weighted mean slope
    if sum_weights == 0:
        return 0.
    else:
        return sum_slopes / sum_weights


def yint(ramp, read_noise, rejection_threshold):
    """
    Find CRs/jumps in a pixel ramp using the y-intercept method.
    """

    good_ends = []
    new_CR = True
    iteration = 0

    # Iterate until no new CRs have been detected
    while new_CR:

        new_CR = False
        iteration += 1

        # Loop over the semi-ramps, looking for outliers in each
        for semiramp in range(ramp.num_semiramps):

            # Get the start and end positions for this semi-ramp
            start = ramp.nodes[semiramp][0]
            end = ramp.nodes[semiramp][1]
            groups = end - start

            # Check to see if this semi-ramp is too short to work with
            if (groups) < 3: good_ends.append(end)

            # If this semi-ramp has already checked out good, then skip it
            if end in good_ends: continue

            # Copy the times and counts for this semi-ramp to working arrays
            times = ramp.times[start:end]
            counts = ramp.counts[start:end]

            # Compute the slopes and y-intercepts for all intervals
            # in this semi-ramp
            slopes, slope_errs, yints, yint_errs = \
                fit_semiramp(times, counts, read_noise)

            # Compute the weighted average of the slope pairs for each
            # interval. Note that the only reason for doing this is to have
            # a signal level from which to estimate the Poisson noise.
            avg_slopes = average_slope_pairs(slopes, slope_errs)

            # Compute photon noise for each average slope value
            pnoise = np.sqrt(np.abs(avg_slopes) * ramp.group_time)

            # Compute the expected uncertainty for each interval based
            # on the combined Poisson and read noise
            yerr_exp = np.empty(len(pnoise))
            for i in range(len(pnoise)):

                # Photon noise and y-intercept uncertainties added
                # in quadrature
                yerr_exp[i] = sqrt(pnoise[i] * pnoise[i] +
                                 yint_errs[i][0] * yint_errs[i][0] +
                                 yint_errs[i][1] * yint_errs[i][1])

            # Compute differences in adjacent y-intercepts for each interval
            ydiff = np.abs(np.diff(yints))

            # Scale the differences by the expected uncertainties
            ratio = ydiff[:, 0] / yerr_exp

            # Check for an outlier that is above rejection threshold
            candidate = ratio.argmax()
            if ratio[candidate] > rejection_threshold:

                # Save the position of the outlier and break
                # this semi-ramp into two pieces
                new_CR = True
                group = start + candidate + 1 # absolute group number in whole ramp
                ramp.groups_with_CR.append(group)
                ramp.nodes[semiramp][1] = group
                if semiramp + 1 < ramp.num_semiramps:
                    ramp.nodes = np.insert(ramp.nodes, semiramp + 1,
                                            [group, start + groups], axis=0)
                else:
                    ramp.nodes = np.append(ramp.nodes,
                                            [[group, start + groups]], axis=0)
                ramp.num_semiramps += 1
                ramp.dqs[group] = ramp.dqs[group] | dqflags.group['JUMP_DET']

                # Break out of loop over semi-ramps and start over
                break

            else:

                # Save the position of this good semi-ramp
                good_ends.append(end)

    return


def fit_semiramp(times, counts, readnoise):
    """
    Fit the slopes and y-intercepts for all possible sample
    combinations within an input semi-ramp
    """

    # Create empty arrays for the results
    groups = len(times)
    slopes = np.empty((groups - 1, 2))
    slope_errs = np.empty((groups - 1, 2))
    yints = np.empty((groups - 1, 2))
    yint_errs = np.empty((groups - 1, 2))

    # Compute fits for the first interval in this semi-ramp.
    #
    # The results for the left-hand side of the first sample pair
    # are set so that the y-intercept is just the counts in the
    # first sample.
    (slopes[0][0], slope_errs[0][0], yints[0][0], yint_errs[0][0]) = \
            0., 999999., counts[0], readnoise

    # The fit for the right-hand side is computed with the time array
    # shifted so that the first point has a positive time value equal
    # to one group time, thus effectively extrapolating the y-intercept
    # backwards one group to the time of the previous sample, matching
    # what was done above for the left-hand side of the pair.
    (slopes[0][1], slope_errs[0][1], yints[0][1], yint_errs[0][1]) = \
            fit_line(times[1:] - times[0], counts[1:],
                     readnoise, random=True)

    # Compute fits for last interval in this semi-ramp.
    #
    # The fit for the left-hand side of the sample pair is computed with
    # the time array shifted so that the last point has a negative value
    # equal to one group time, thus effectively extrapolating the
    # y-intercept forward one group to the time of the next sample,
    # matching what will be done below for the right-hand side of the pair.
    (slopes[groups - 2][0], slope_errs[groups - 2][0],
            yints[groups - 2][0], yint_errs[groups - 2][0]) = \
            fit_line(times[:groups - 1] - times[groups - 1], counts[:groups - 1],
                     readnoise, random=True)

    # The results for the right-hand side of the last sample pair are
    # set so that the y-intercept is just the counts in the last sample.
    (slopes[groups - 2][1], slope_errs[groups - 2][1],
            yints[groups - 2][1], yint_errs[groups - 2][1]) = \
            0., 999999., counts[groups - 1], readnoise

    # Compute fits for all intermediate intervals
    for group in range(1, groups - 2):

        # For the left-hand side of each sample pair, shift the time array
        # so that the last point has a negative value equal to one group
        # time, thus effectively extrapolating the y-intercept forward one
        # group to the time of the next sample, matching what will be done
        # below for the right-hand side of the pair.
        (slopes[group][0], slope_errs[group][0],
            yints[group][0], yint_errs[group][0]) = \
            fit_line(times[:group + 1] - times[group + 1], counts[:group + 1],
                     readnoise, random=True)
            # Use the following line to shift the y-intercept to the
            # left side of the pair.
            #fit_line(times[:group+1]-times[group],counts[:group+1], random=True)

        # For the right-hand side of the sample pair, shift the time array
        # so that the first point has a time of zero, thus setting the
        # y-intercept to correspond to that sample.
        (slopes[group][1], slope_errs[group][1],
            yints[group][1], yint_errs[group][1]) = \
            fit_line(times[group + 1:] - times[group + 1], counts[group + 1:],
                     readnoise, random=True)
            # Use the following line to shift the y-intercept to the
            # left side of the pair.
            #fit_line(times[group+1:]-times[group],counts[group+1:], random=True)

    return slopes, slope_errs, yints, yint_errs


def average_slope_pairs(slopes, slope_errs):
    """
    Compute the weighted average of the slopes for each sample pair
    within a semi-ramp.
    """

    weight = 1. / (slope_errs * slope_errs)
    avg_slopes = np.sum(slopes * weight, axis=1) / np.sum(weight, axis=1)

    return avg_slopes


def fit_line(xx, yy, readnoise, random=False):
    """
    Fit a 1st-order polynomial to a set of input x/y ramp values.
    The fitting is performed using matrix arithmetic, taking
    into account the full co-variance matrix of the readouts.
    """

    # Initial slope estimate, from which Poisson noise is estimated
    slope, yint = np.polyfit(xx, yy, 1)
    photon_noise = sqrt(abs(slope) * (xx[1] - xx[0]))
    pn2 = photon_noise * photon_noise

    nn = len(xx)
    if nn < 3:
        mm = slope
        bb = yint
        mm_err = sqrt(pn2 / 2.)
        bb_err = readnoise
        return mm, mm_err, bb, bb_err

    # Initialize matrices
    YY = np.asmatrix(yy).T
    AA = np.asmatrix([np.ones(nn), xx]).T
    guts = [np.zeros(nn)] * nn
    PP = np.asmatrix(guts)

    # Photon noise (correlated) (symmetric matrix)
    for j in range(nn):
        for i in range(j, nn):
            PP[i, j] = PP[j, i] = pn2 * (j + 1)

    # Read noise (uncorrelated) (on the diagonal)
    RR = np.asmatrix(np.identity(nn) * readnoise * readnoise)

    # Full covariance matrix
    CC = PP + RR

    # Compute intercept and slope by solving matrices
    part1 = AA.T * (CC.I * AA)
    part2 = AA.T * (CC.I * YY)
    XX = part1.I * part2
    bb, mm = XX[0, 0], XX[1, 0]

    # Calculate uncertainties on slope and intercept
    mm_err = sqrt(part1.I[1, 1])
    if random:
        part1 = AA.T * (RR.I * AA)
    bb_err = sqrt(part1.I[0, 0])

    return mm, mm_err, bb, bb_err
