import numpy as np
import sys
import random as rand
import subprocess

random = rand.SystemRandom()

weighers = ["UniformWeigher",
            "InversePredicateFrequencyWeigher",
            "PredicateFrequencyWeigher",
            "ObjectFrequencyWeigher",
            "InverseObjectFrequencyWeigher",
            "PredicateObjectFrequencyWeigher",
            "InversePredicateObjectFrequencyWeigher",
            "InverseObjectFrequencyWeigherSplitDown",
            "PushDownWeigherMap"]

allCombs = np.array(np.meshgrid(weighers, weighers)).T.reshape(-1,2)
alphas = np.arange(0.4, 1.0, 0.05)
epsilons = np.power(np.array([10.]), np.array([-6,-5,-4,-3]))
normalize = [True,False]
onlyEntities = [True,False]

def GenCall():
    call = ["/usr/bin/time", "-v","./main", "-l"]
    n = random.choice(normalize)
    o = random.choice(onlyEntities)
    alpha = random.choice(alphas)
    eps = random.choice(epsilons)
    ws = random.choice(allCombs)
    if n == True:
        call.append("-n")
    if o == True:
        call.append("--only-entities")
    call.append("-a")
    call.append(str(alpha))
    call.append("-e")
    call.append(str(eps))
    call.append("-f")
    call.append(ws[0])
    call.append("-b")
    call.append(ws[1])
    return call



if __name__ == '__main__':
    alreadyCalled = [["",""]]
    for sample in range(60):
        call = ["",""]
        while call in alreadyCalled:
            call = GenCall()
        alreadyCalled.append(call)
        print(call)
        sys.stdout.flush()
        res = subprocess.call(call)
        if res == 0:
            print("succeeded\n")
        else:
            print("failed\n")
        sys.stdout.flush()
        sys.stderr.flush()
