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

def get_detectors(filename):
    tree = et.parse(filename)
    proposal = tree.getroot()
    ObservationNum = 0
    detectors = []
    ProposalNum = filename.split('.')[0]
    #
    # Force ProposalNum to be eactly 5 characters long by either truncating or padding
    if len(ProposalNum) > 5:
        ProposalNum = ProposalNum[-5:]
    elif len(ProposalNum) < 5:
        ProposalNum = ProposalNum.zfill(5)
    for observation in proposal.getchildren():
        ObservationNum = ObservationNum + 1
        VisitNum = 0
        for visit in observation.getchildren():
            VisitNum = VisitNum + 1
            VisitGroupNum = 0
            for visitgroup in visit.getchildren():
                VisitGroupNum = VisitGroupNum + 1
                ParallelSequenceIDNum = 0
                for parallelsequenceid in visitgroup.getchildren():
                    ParallelSequenceIDNum = ParallelSequenceIDNum + 1
                    ActivityNum = 0
                    for activity in parallelsequenceid.getchildren():
                        ActivityNum = ActivityNum + 1
                        ExposureNum = 0
                        for exposure in activity.getchildren():
                            ExposureNum = ExposureNum + 1
                            DetectorNum = 0
                            for detector in exposure.getchildren():
                                subarray = detector.find('subarray').text
                                base = detector.find('base').text
                                exp_type = detector.find('exp_type').text
                                id = ProposalNum + '%03d%03d%02d%d%02d%05d' % (ObservationNum,
                                                                               VisitNum,
                                                                               VisitGroupNum,
                                                                               ParallelSequenceIDNum,
                                                                               ActivityNum,
                                                                               ExposureNum)
                                detectors.append({'base': base,
                                                  'subarray': subarray,
                                                  'exp_type': exp_type,
                                                  'id': id})
    return detectors
