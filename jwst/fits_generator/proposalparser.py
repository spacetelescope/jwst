#
#  Proposal parser.  A proposal follows this structure:
#
#  <Proposal>
#    <Observation>
#      <Visit>
#        <VisitGroup>
#          <ParallelSequenceID>
#            <Activity>
#              <Exposure>
#                <Detector>
#                  <base></base>
#                  <subarray></subarray>
#                  <exp_type></exp_type>
#                </Detector>
#              </Exposure>
#            </Activity>
#          </ParallelSequenceID>
#        </VisitGroup>
#      </Visit>
#     </Observation>
#   </Proposal>
#
# Each element can occur 1 or more times.
#

import xml.etree.ElementTree as et

class XMLProposal(object):
    def __init__(self, filename):
        self.tree = et.parse(filename)
        self.proposal = self.tree.getroot()
        self.proposal_id = filename.split('.')[0]

    def get_proposal_number(self):
        pass

    def get_observations(self):
        self.observations = self.proposal.getchildren()
        n = 0
        for observation in self.observations:
            n = n + 1
            observation.value = '%03d' % n
            observation.id = ''.join((self.proposal_id, '_', observation.value))
        
        return self.observations

    def get_visits(self):
        self.get_observations()
        self.visits = []
        for observation in self.observations:
            visits = observation.getchildren()
            n = 0
            for visit in visits:
                n = n + 1
                visit.value = '%03d' % n
                visit.id = ''.join((observation.id, '_', visit.value))
                self.visits.append(visit)

        return self.visits

    def get_visit_groups(self):
        self.get_visits()
        self.visit_groups = []
        for visit in self.visits:
            visit_groups = visit.getchildren()
            n = 0
            for visit_group in visit_groups:
                n = n + 1
                visit_group.value = '%02d' % n
                visit_group.id = ''.join((visit.id, '_', visit_group.value))
                self.visit_groups.append(visit_group)

        return self.visit_groups

    def get_parallel_sequence_ids(self):
        self.get_visit_groups()
        self.parallel_sequence_ids = []
        for visit_group in self.visit_groups:
            parallel_sequence_ids = visit_group.getchildren()
            n = 0
            for parallel_sequence_id in parallel_sequence_ids:
                n = n + 1
                parallel_sequence_id.value = '%d' % n
                parallel_sequence_id.id = ''.join((visit_group.id, parallel_sequence_id.value))
                self.parallel_sequence_ids.append(parallel_sequence_id)

        return self.parallel_sequence_ids

    def get_activities(self):
        self.get_parallel_sequence_ids()
        self.activities = []
        for parallel_sequence_id in self.parallel_sequence_ids:
            activities = parallel_sequence_id.getchildren()
            n = 0
            for activity in activities:
                n = n + 1
                activity.value = '%02d' % n
                activity.id = ''.join((parallel_sequence_id.id, activity.value))
                self.activities.append(activity)

        return self.activities

    def get_exposures(self):
        self.get_activities()
        self.exposures = []
        for activity in self.activities:
            exposures = activity.getchildren()
            n = 0
            for exposure in exposures:
                n = n + 1
                exposure.value = '%05d' % n
                exposure.id = ''.join((activity.id, '_', exposure.value))
                self.exposures.append(exposure)
                
        return self.exposures

    def get_detectors(self):
        self.get_exposures()
        self.detectors = []
        for exposure in self.exposures:
            detectors = exposure.getchildren()
            n = 0
            for detector in detectors:
                n = n + 1
                detector.subarray = detector.find('subarray').text
                detector.base = detector.find('base').text
                detector.exp_type = detector.find('exp_type').text
                detector.id = exposure.id
                self.detectors.append(detector)
                
        return self.detectors
