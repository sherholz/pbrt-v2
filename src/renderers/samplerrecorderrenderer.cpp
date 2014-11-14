
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


// renderers/samplerrecorderrenderer.cpp*
#include "stdafx.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include "renderers/samplerrecorderrenderer.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "integrators/path.h"


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

// SamplerRecorderRendererTask Definitions
void SamplerRecorderRendererTask::Run() {
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
    RayDifferential ray;
    Spectrum T;
	vector<Spectrum> L(3,0.);
    Intersection isect;

    // Get samples from _Sampler_ and update image
	int itr_count = 0;
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {

		if(sampleIndex>=sampleCount) return;
		itr_count++;
        // Generate camera rays and compute radiance along rays
       
        // Find camera ray for _sample[sampleIndex]_
        PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[sampleIndex]);
        float rayWeight = cameras.at(0)->GenerateRayDifferential(samples[sampleIndex], &ray);
        ray.ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
        PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[sampleIndex], &rays[sampleIndex], rayWeight);

        // Evaluate radiance along camera ray
        for(int j=0; j<L.size();j++) PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[sampleIndex], &samples[sampleIndex]);
        if (visualizeObjectIds) {
            if (rayWeight > 0.f && scene->Intersect(ray, &isect)) {
                // random shading based on shape id...
                uint32_t ids[2] = { isect.shapeId, isect.primitiveId };
                uint32_t h = hash((char *)ids, sizeof(ids));
                float rgb[3] = { float(h & 0xff), float((h >> 8) & 0xff),
                                    float((h >> 16) & 0xff) };
                for(int j=0; j<L.size();j++)L.at(j) = Spectrum::FromRGB(rgb);
                for(int j=0; j<L.size();j++)L.at(j) /= 255.f;
            }
            else
                for(int j=0; j<L.size();j++)L.at(j) = 0.f;
        }
        else {
			if (rayWeight > 0.f){
				L = renderer->Li_separated(scene, ray, &samples[sampleIndex], rng, arena, &isect, &T);
				for(int j=0; j<L.size();j++) L.at(j)*=rayWeight;
			} else {
				for(int j=0; j<L.size();j++)L.at(j) = 0.f;
				T = 1.f;
			}

			for(int j=0; j<L.size();j++)
			{
				// Issue warning if unexpected radiance value returned
				if (L.at(j).HasNaNs()) {
					Error("Not-a-number radiance value returned "
							"for image sample.  Setting to black.");
					L.at(j) = Spectrum(0.f);
				}
				else if (L.at(j).y() < -1e-5) {
					Error("Negative luminance value, %f, returned "
							"for image sample.  Setting to black.", L.at(j).y());
					L.at(j) = Spectrum(0.f);
				}
				else if (isinf(L.at(j).y())) {
					Error("Infinite luminance value returned "
							"for image sample.  Setting to black.");
					L.at(j) = Spectrum(0.f);
				}
			}
        }
        for(int j=0; j<L.size();j++) PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&ray, &samples[sampleIndex], &L.at(j));

		for(int j=0; j<L.size();j++)
		{
			// Report sample results to _Sampler_, add contributions to image
			if (sampler->ReportResults(&samples[sampleIndex], &ray, &L.at(j), &isect, 1))
			{ 
				PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[sampleIndex], &ray, &L.at(j),&T);

				// use SetSample instead of AddSample to get only the current sample contribution
				cameras.at(j)->film->SetSample(samples[sampleIndex], L.at(j));
				PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
			}
		}
        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }

    // Clean up after _SamplerRendererTask_ is done with its image region
    for(int j=0; j<L.size();j++) cameras.at(j)->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    delete[] samples;
    reporter.Update();

    PBRT_FINISHED_RENDERTASK(taskNum);
}

void SamplerAccumulationRendererTask::Run() {
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
    vector<Spectrum *>Ls(3,0);
	for(int i=0; i<Ls.size();i++)Ls.at(i) = new Spectrum[maxSamples];

	Spectrum * Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = cameras.at(0)->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

            // Evaluate radiance along camera ray
            for(int j=0; j<Ls.size();j++) PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
            if (visualizeObjectIds) {
                if (rayWeight > 0.f && scene->Intersect(rays[i], &isects[i])) {
                    // random shading based on shape id...
                    uint32_t ids[2] = { isects[i].shapeId, isects[i].primitiveId };
                    uint32_t h = hash((char *)ids, sizeof(ids));
                    float rgb[3] = { float(h & 0xff), float((h >> 8) & 0xff),
                                     float((h >> 16) & 0xff) };
                    for(int j=0; j<Ls.size();j++)Ls.at(j)[i] = Spectrum::FromRGB(rgb);
                    for(int j=0; j<Ls.size();j++)Ls.at(j)[i] /= 255.f;
                }
                else
                    for(int j=0; j<Ls.size();j++)Ls.at(j)[i] = 0.f;
            }
            else {
				if (rayWeight > 0.f){
					vector<Spectrum> tmp = renderer->Li_separated(scene, rays[i], &samples[i], rng, arena, &isects[i], &Ts[i]);
					for(int j=0; j<Ls.size();j++)Ls.at(j)[i] = rayWeight * tmp.at(j);
				}else {
					for(int j=0; j<Ls.size();j++)Ls.at(j)[i] = 0.f;
					Ts[i] = 1.f;
				}

				// Issue warning if unexpected radiance value returned
				for(int j=0; j<Ls.size();j++)
				{
					if (Ls.at(j)[i].HasNaNs()) {
						Error("Not-a-number radiance value returned "
							  "for image sample.  Setting to black.");
						Ls.at(j)[i] = Spectrum(0.f);
					}
					else if (Ls.at(j)[i].y() < -1e-5) {
						Error("Negative luminance value, %f, returned "
							  "for image sample.  Setting to black.", Ls.at(j)[i].y());
						Ls.at(j)[i] = Spectrum(0.f);
					}
					else if (isinf(Ls.at(j)[i].y())) {
						Error("Infinite luminance value returned "
							  "for image sample.  Setting to black.");
						Ls.at(j)[i] = Spectrum(0.f);
					}
				}
            }
            for(int j=0; j<Ls.size();j++) PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls.at(j)[i]);
        }

		for(int j=0; j<Ls.size();j++)
		{
			// Report sample results to _Sampler_, add contributions to image
			if (sampler->ReportResults(samples, rays, Ls.at(j), isects, sampleCount))
			{
				for (int i = 0; i < sampleCount; ++i)
				{
					if(j==0)PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls.at(j)[i], &Ts[i]);
					cameras.at(j)->film->AddSample(samples[i], Ls.at(j)[i]);
					if(j==0)PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
				}
			}
		}

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }

    // Clean up after _SamplerLightSeparationRendererTask_ is done with its image region
    for(int j=0; j<Ls.size();j++) cameras.at(j)->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);

    delete sampler;
    delete[] samples;
    delete[] rays;
    for(int j=0; j<Ls.size();j++) delete[] Ls.at(j);
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}

// SamplerRecorderRenderer Method Definitions
SamplerRecorderRenderer::SamplerRecorderRenderer(Sampler *s, vector<Camera *>c,
                                 PathIntegrator *pi, VolumeIntegrator *vi,
                                 bool visIds) {
    sampler = s;
    cameras = c;
	pathIntegrator = pi;
    volumeIntegrator = vi;
    visualizeObjectIds = visIds;
}


SamplerRecorderRenderer::~SamplerRecorderRenderer() {
    delete sampler;
    for(int j=0; j<cameras.size();j++) delete cameras.at(j);
    delete pathIntegrator;
    delete volumeIntegrator;
}


void SamplerRecorderRenderer::Render(const Scene *scene) {
   
	int sampleCount = sampler->MaximumSampleCount();
	PBRT_FINISHED_PARSING();

	// Allow integrators to do preprocessing for the scene
	PBRT_STARTED_PREPROCESSING();
	for(int j=0; j<cameras.size();j++)pathIntegrator->Preprocess(scene, cameras.at(j), this);
	for(int j=0; j<cameras.size();j++)volumeIntegrator->Preprocess(scene, cameras.at(j), this);
	PBRT_FINISHED_PREPROCESSING();

	PBRT_STARTED_RENDERING();
	// Allocate and initialize _sample_
	Sample *sample = new Sample(sampler, pathIntegrator,
								volumeIntegrator, scene);

	// Compute number of _SamplerRendererTask_s to create for rendering
	int nPixels = cameras.at(0)->film->xResolution * cameras.at(0)->film->yResolution;
	int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
	nTasks = RoundUpPow2(nTasks);

	// write out an image after each sample
	for(int s = 0; s<sampleCount;s++)
	{
		std::stringstream ss;
		ss<<"Rendering Sample "<<s+1;
		ProgressReporter reporter(nTasks, ss.str());

		vector<Task *> renderTasks;

		for (int i = 0; i < nTasks; ++i)
			renderTasks.push_back(new SamplerRecorderRendererTask(scene, this, cameras,
														  reporter, sampler, sample, 
														  visualizeObjectIds, 
														  nTasks-1-i, nTasks,s));
		EnqueueTasks(renderTasks);
		WaitForAllTasks();
		for (uint32_t i = 0; i < renderTasks.size(); ++i)
			delete renderTasks[i];

		// write out samples for direct light
		ss.str("");
		ss<<"direct_" <<std::setfill('0') << std::setw(3)<<s+1;
		cameras.at(0)->film->WriteImageWithSuffix(ss.str());

		// write out samples for indirect light
		ss.str("");
		ss<<"indirect_"<<std::setfill('0') << std::setw(3)  <<s+1;
		cameras.at(1)->film->WriteImageWithSuffix(ss.str());

		// write out samples for direct light
		ss.str("");
		ss <<"combined_"<<std::setfill('0') << std::setw(3) <<s+1;
		cameras.at(2)->film->WriteImageWithSuffix(ss.str());

		reporter.Done();
	}

	/*
	// write accumulated image
	ProgressReporter reporter(nTasks, "Rendering accumulated Samples");
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i)
        renderTasks.push_back(new SamplerAccumulationRendererTask(scene, this, cameras,
                                                      reporter, sampler, sample, 
                                                      visualizeObjectIds, 
                                                      nTasks-1-i, nTasks));
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
        delete renderTasks[i];
    reporter.Done();

	// write out accumulated image
	cameras.at(0)->film->WriteImageWithSuffix("accum_direct");
	cameras.at(1)->film->WriteImageWithSuffix("accum_indirect");
	cameras.at(2)->film->WriteImageWithSuffix("accum_combined");
	*/
	// Clean up after rendering and store final image
	PBRT_FINISHED_RENDERING();
	delete sample;

}

vector<Spectrum> SamplerRecorderRenderer::Li_separated(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());

    // Allocate local variable for _isect_ if needed

    Intersection localIsect;
    if (!isect) isect = &localIsect;
    vector<Spectrum> Li (3,0.);
    if (scene->Intersect(ray, isect))
		Li = pathIntegrator->Li_separate(scene, this, ray, *isect, sample,
                                   rng, arena);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
			for (uint32_t j = 0; j < Li.size(); ++j)
				Li.at(j) += scene->lights[i]->Le(ray);
    }

    return Li;
}

Spectrum SamplerRecorderRenderer::Li(const Scene *scene,
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
    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
                                        T, arena);
    return *T * Li + Lvi;
}


Spectrum SamplerRecorderRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}


