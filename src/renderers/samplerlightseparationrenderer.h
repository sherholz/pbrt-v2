
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_LIGHT_SEPARATION_RENDERER_H
#define PBRT_RENDERERS_LIGHT_SEPARATION_RENDERER_H

// renderers/samplerlightseparationrenderer.h*
#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"

class PathIntegrator;

enum IlluminationType
{
	DIRECT_ILLUMINATION,
	INDIRECT_ILLUMINATION
};

// SamplerLightSeparationRenderer Declarations
class SamplerLightSeparationRenderer : public Renderer {
public:
    // SamplerLightSeparationRenderer Public Methods
    SamplerLightSeparationRenderer(Sampler *s, vector<Camera *>c, PathIntegrator *pi,
                    VolumeIntegrator *vi, bool visIds);
    ~SamplerLightSeparationRenderer();
    void Render(const Scene *scene);

	// returns direct (index 0), indirect (index 1) and combined (index 2) illumination 
    vector<Spectrum> Li_separated(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL) const;

	Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL) const;
	Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    
private:
    // SamplerLightSeparationRenderer Private Data
    bool visualizeObjectIds;
    Sampler *sampler;
    vector<Camera *>cameras;
    PathIntegrator *pathIntegrator;
    VolumeIntegrator *volumeIntegrator;
};



// SamplerLightSeparationRendererTask Declarations
class SamplerLightSeparationRendererTask : public Task {
public:
    // SamplerLightSeparationRendererTask Public Methods
    SamplerLightSeparationRendererTask(const Scene *sc, SamplerLightSeparationRenderer *ren, vector<Camera *>c,
                        ProgressReporter &pr, Sampler *ms, Sample *sam, 
                        bool visIds, int tn, int tc, IlluminationType it)
      : reporter(pr)
    {
        scene = sc; renderer = ren; cameras = c; mainSampler = ms;
        origSample = sam; visualizeObjectIds = visIds; taskNum = tn; taskCount = tc;
		illuminationType = it;
    }
    void Run();
private:
    // SamplerLightSeparationRendererTask Private Data
    const Scene *scene;
    const SamplerLightSeparationRenderer *renderer;
    vector<Camera *>cameras;
    Sampler *mainSampler;
    ProgressReporter &reporter;
    Sample *origSample;
    bool visualizeObjectIds;
    int taskNum, taskCount;
	IlluminationType illuminationType;
};



#endif // PBRT_RENDERERS_SAMPLER_LIGHT_SEPARATION_RENDERER_H
