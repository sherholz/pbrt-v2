
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


// renderers/samplerrenderer.cpp*
#include "stdafx.h"
#include "renderers/samplerrenderer.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "integrators/path.h"
#include "integrators/single.h"

static uint32_t hash(char *key, uint32_t len)
{
    uint32_t hash = 0, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
} 

// SamplerRendererTask Definitions
void SamplerRendererTask::Run() {
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _SamplerRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }

    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

            // Evaluate radiance along camera ray
            PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
            
            if (rayWeight > 0.f)
                Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
                                                 arena, &isects[i], &Ts[i]);
            else {
                Ls[i] = 0.f;
                Ts[i] = 1.f;
            }

            // Issue warning if unexpected radiance value returned
            if (Ls[i].HasNaNs()) {
                Error("Not-a-number radiance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            else if (Ls[i].y() < -1e-5) {
                Error("Negative luminance value, %f, returned "
                      "for image sample.  Setting to black.", Ls[i].y());
                Ls[i] = Spectrum(0.f);
            }
            else if (isinf(Ls[i].y())) {
                Error("Infinite luminance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            
            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }

        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
                PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                camera->film->AddSample(samples[i], Ls[i]);
                PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
            }
        }

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }

    // Clean up after _SamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}


	
void SamplerRendererTaskTmp::Run() {
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _SamplerLightSeparationRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }

    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
	
	SpectrumContainer* Ls = new SpectrumContainer[maxSamples];
	for(int i=0; i<maxSamples;i++) Ls[i].SetRenderPassTypes(camera->GetRenderPassTypes());

	Spectrum * Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
			// Clear the SpectrumContainer
			Ls[i].Clear(Spectrum(0.));

            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

            // Evaluate radiance along camera ray
            PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
			
			bool intersected = scene->Intersect(rays[i], &isects[i]);

            if (Ls[i].ContainsSpectrum(PRIMITIVES)) {
                if (rayWeight > 0.f && intersected) {
                    // random shading based on shape id...
					uint32_t ids[2] = { isects[i].shapeId, isects[i].primitiveId };
                    uint32_t h = hash((char *)ids, sizeof(ids));
                    float rgb[3] = { float(h & 0xff), float((h >> 8) & 0xff), float((h >> 16) & 0xff) };
                    Ls[i][PRIMITIVES] = Spectrum::FromRGB(rgb);
                    Ls[i][PRIMITIVES] = Ls[i][PRIMITIVES]/255.f;
                }
                else{
                    Ls[i][PRIMITIVES] = Spectrum(0.);
				}
            }
           
			if (rayWeight > 0.f){
				renderer->Li_separate(scene, rays[i], &samples[i], rng, arena, Ls[i], &isects[i], intersected, &Ts[i]);

				// only scale radiance value by ray weight
				Ls[i][DIRECT]*=rayWeight;
				Ls[i][INDIRECT]*=rayWeight;
				Ls[i][BEAUTY]*=rayWeight;
			}else {
				Ls[i] = Spectrum(0.);
				Ts[i] = 1.f;
			}

			// Issue warning if unexpected radiance value returned
			Spectrum spec = Ls[i].GetFirstSpectrum();
			while(Ls[i].HasNextSpectrum())
			{					
				if (spec.HasNaNs()) {
					Error("Not-a-number value returned for image sample. Setting to black.");
					spec = Spectrum(0.f);
				}
				else if (spec.y() < -1e-5) {
					Error("Negative value, %f, returned for image sample. Setting to black.", spec.y());
					spec = Spectrum(0.f);
				}
				else if (isinf(spec.y())) {
					Error("Infinite value returned for image sample. Setting to black.");
					spec = Spectrum(0.f);
				}

				spec = Ls[i].GetNextSpectrum();
			}
        }

        PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i].GetFirstSpectrum());

		for(size_t i=0; i<camera->GetNumRenderPasses();i++)
		{
			Film* activePass = camera->GetRenderPass(i);
			RenderPassType activePassType = activePass->GetRenderPassType();

			// Report sample results to _Sampler_, add contributions to image
			if (sampler->ReportResults(samples, rays, &Ls[0][activePassType], isects, sampleCount))
			{
				for (int s = 0; s < sampleCount; ++s)
				{
					PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[s], &rays[s], &Ls[s][activePassType], &Ts[s]);
					activePass->AddSample(samples[s], Ls[s][activePassType]);
					PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
				}
			}
			activePass->UpdateDisplay(sampler->xPixelStart, sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
		}

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }

    // Clean up after _SamplerLightSeparationRendererTask_ is done with its image region
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ts;
	delete[] Ls;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}


// SamplerRenderer Method Definitions
SamplerRenderer::SamplerRenderer(Sampler *s, Camera *c,PathIntegrator *pi, 
								 SingleScatteringIntegrator *ssi) {
    sampler = s;
    camera = c;
    pathIntegrator = pi;
    singleScatteringIntegrator = ssi;
}


SamplerRenderer::~SamplerRenderer() {
    delete sampler;
    delete camera;
    delete pathIntegrator;
    delete singleScatteringIntegrator;
}


void SamplerRenderer::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    pathIntegrator->Preprocess(scene, camera, this);
	singleScatteringIntegrator->Preprocess(scene, camera, this);
    PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, pathIntegrator,
								singleScatteringIntegrator, scene);

    // Compute number of _SamplerRendererTask_s to create for rendering
    int nPixels = camera->film->xResolution * camera->film->yResolution;
    int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
    nTasks = RoundUpPow2(nTasks);
    ProgressReporter reporter(nTasks, "Rendering");
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i)
        renderTasks.push_back(new SamplerRendererTaskTmp(scene, this, camera,
                                                      reporter, sampler, sample, 
                                                      nTasks-1-i, nTasks));
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
        delete renderTasks[i];
    reporter.Done();
    PBRT_FINISHED_RENDERING();

	// write out all render passes
	for(size_t i=0; i<camera->GetNumRenderPasses();i++){
		camera->GetRenderPass(i)->WriteImage();
	}

    // Clean up after rendering and store final image
    delete sample;
}


void SamplerRenderer::Li_separate(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
		MemoryArena &arena, SpectrumContainer& L_io, Intersection *isect, bool intersected, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());

    // Allocate local variable for _isect_ if needed
	Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    
    if (intersected)//scene->Intersect(ray, isect))
		// result is stored in L_io
		pathIntegrator->Li_separate(scene, this, ray, *isect, sample, rng, arena, L_io);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i){
			Spectrum Le = scene->lights[i]->Le(ray);
			L_io[BEAUTY] += Le;
			L_io[DIRECT] += Le;
			L_io[INDIRECT] += Le;
		}
    }
	
	SpectrumContainer Lv_io(camera->GetRenderPassTypes());
	singleScatteringIntegrator->Li_separate(scene, this, ray, sample, rng,T, arena, Lv_io);

	L_io[BEAUTY] += Lv_io[BEAUTY];
	L_io[DIRECT] += Lv_io[DIRECT];
	L_io[INDIRECT] += Lv_io[INDIRECT];
}

Spectrum SamplerRenderer::Li(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect))
        Li = pathIntegrator->Li(scene, this, ray, *isect, sample,
                                   rng, arena);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           Li += scene->lights[i]->Le(ray);
    }
    Spectrum Lvi = singleScatteringIntegrator->Li(scene, this, ray, sample, rng,
                                        T, arena);
    return *T * Li + Lvi;
}


Spectrum SamplerRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return singleScatteringIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}


